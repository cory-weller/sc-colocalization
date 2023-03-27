#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

DISEASE_GWAS_FILES <- list.files(path='/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats/',
                                    pattern='*.formatted.tsv',
                                    full.names=T)

RECOMBINATION_BED_FILE <- '/gpfs/gsfs9/users/wellerca/sc-colocalization/data/genetic_map_hg38_withX.txt.gz'

GWAS_P_THRESHOLD <- 5e-8
cM_THRESHOLD <- 0.2

import_recombination_bed <- function(filename) {
    # Read in .bed file containing recombination rates
    # Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
    bed <- fread(file = RECOMBINATION_BED_FILE,
                header = TRUE,
                showProgress = FALSE,
                col.names = c("CHR","BP","c","cM")
    )
    bed[, c := NULL]
    setkey(bed, CHR, BP)
    bed_extension <- bed[, list('BP'=max(BP)*1.1, 'cM'=max(cM)), by=CHR]    # extend ends of chr by 10%, with same cM
    bed <- unique(rbindlist(list(bed, bed_extension)))
    setkey(bed, CHR, BP)
    return(bed)
}


build_pos_to_cM <- function(recombination_bed_table) {
# Generate functions to translate base pair position to units of cM
    #cM_to_POS <- new.env()
    env1 <- new.env()
    for(chr in unique(bed[,CHR])) {
        #cM_to_POS[[as.character(chr)]] <- approxfun(y=c(0, bed[CHR==chr][,BP]), x=c(0,bed[CHR==chr][,cM]), ties='ordered')
        env1[[as.character(chr)]] <- approxfun(x=c(0, bed[CHR==chr][,BP]), y=c(0,bed[CHR==chr][,cM]), ties='ordered')
    }
    return(env1)
}


build_cM_to_pos <- function(recombination_bed_table) {
# Generate function to translate cM to base pair position
    env1 <- new.env()
    for(chr in unique(bed[,CHR])) {
        env1[[as.character(chr)]] <- approxfun(y=c(0, bed[CHR==chr][,BP]), x=c(0,bed[CHR==chr][,cM]), ties='ordered')
    }
    return(env1)
}



get_gwas_filestem <- function(filename) {
    return(strsplit(basename(filename), split='\\.')[[1]][1])
}


get_gwas_hits <- function(gwas_filename) {
    dat <- fread(gwas_filename)
    setnames(dat, toupper(colnames(dat)))
    dat <- dat[P <= GWAS_P_THRESHOLD]
    setkey(dat, CHR, BP)


    desired_cols <- c('SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P')
    dat <- dat[, .SD, .SDcols=desired_cols]

    o <- foreach(chr=1:21, .combine='rbind') %do% {
        dat.sub <- dat[CHR==chr]
        dat.sub[, cM := pos_to_cM[[as.character(chr)]](BP)]
        dat.sub[, nextSNP := shift(BP, n=-1L), by=CHR]
        dat.sub[, nextSNP_cM := pos_to_cM[[as.character(chr)]](nextSNP)]
        dat.sub[, cM_between_SNPs := nextSNP_cM - cM]
        return(dat.sub)
    }
    return(o)
}

get_clusters <- function(DT, cM_threshold) {
    o <- copy(DT)
    o[, tf:= ifelse(cM_between_SNPs <= cM_threshold, TRUE, FALSE)]    # TRUE if next SNP is within cM_THRESHOLD, otherwise FALSE
    o[ is.na(nextSNP), tf := FALSE]

    # vectors of what cluster start and stop ROWS from o should be
    starts <- unique(c(1, which(o$tf == FALSE)+1))
    stops <- unique(c(which(o$tf == FALSE), nrow(o)))

    o2 <- foreach(start=starts, stop=stops, grp=1:length(starts), .combine='rbind') %do% {
        o.sub <- o[start:stop]
        o.sub[, 'cluster' := grp]
        return(o.sub[])
    }

    o2[, tf := NULL]

    o3 <- o2[, list('bp_start'=min(BP),
              'bp_stop'=max(BP),
              'cM_start'=min(cM),
              'cM_stop'=max(cM),
              'nSNPs_in_cluster'=.N), by=list(CHR,cluster)]
    
    # o3[, cM_width := cM_stop - cM_start]
    # o3[, bp_width := bp_stop - bp_start]
    o3[, 'cM_threshold' := cM_threshold]
    return(o3)

    # # Does not run--code block for manual debugging/wrangling
    # if(FALSE) {
    #     # Assign label for first/last BP of group
    #     o2[o2[, .I[which.min(BP)], by=grp]$V1, ggLabel := BP]
    #     o2[o2[, .I[which.max(BP)], by=grp]$V1, ggLabel := BP]


    #     g <- ggplot(o2, aes(x=cM, y=1, label=ggLabel, color=as.factor(chrGrp))) +
    #     geom_point() +
    #     geom_text_repel() +
    #     facet_grid(CHR~.) +
    #     theme_few()
    # }
    # if(write_output) {
    #     fwrite(o3, file=paste0(gwas_name, '.cluster_windows.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
    # } else {
    #     return(o3)
    # }
}

extend_clusters <- function(DT, cM_to_extend) {
    o <- foreach(chr=1:21, .combine='rbind') %do% {
        dat.sub <- copy(DT[CHR==chr])
        dat.sub[, extended_bp_start := cM_to_pos[[as.character(chr)]](cM_start - 0.5)]
        dat.sub[, extended_bp_stop := cM_to_pos[[as.character(chr)]](cM_start + 0.5)]
        return(dat.sub)
    }
}


bed <- import_recombination_bed(RECOMBINATION_BED_FILE)
pos_to_cM <- build_pos_to_cM(bed)
cM_to_pos <- build_cM_to_pos(bed)


o2 <- foreach(gwas_filename = DISEASE_GWAS_FILES, .combine='rbind') %do% {
    gwas_name <- get_gwas_filestem(gwas_filename)
    signif_snps <- get_gwas_hits(gwas_filename)
    o <- foreach(cM=c(0.01,0.02, 0.03, 0.04, seq(0.05, 1, 0.05)), .combine='rbind') %do% {
        get_clusters(signif_snps, cM)
    }
    o[, 'GWAS' := gwas_name]
}

setkey(o2, GWAS, cM_threshold, CHR, cluster)
desired_col_order <- c('GWAS', 'cM_threshold', 'cluster', 'CHR', 'bp_start', 'bp_stop', 'bp_width',
                    'cM_start', 'cM_stop', 'cM_width','nSNPs_in_cluster')

clusters <- o2[, .SD, .SDcols=desired_col_order]
fwrite(clusters, file='data/clusters.tsv', quote=F, row.names=F, col.names=T, sep='\t')

g <- ggplot(cluster_sizes.long[variable=='N_clusters'] , aes(x=cM_threshold, y=count)) +
    geom_point(shape=21) +
    facet_wrap(~GWAS, nrow=2) +
    theme_few() +
    geom_vline(xintercept=0.2, alpha=0.4, linetype='dashed', color='red') +
    theme(strip.placement = "outside") +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x='cM threshold to cluster adjacent SNPs',
        y='Count',
        title='Number of SNP clusters (across all chromosomes) for colocalization analysis')

ggsave(g, file='data/colocalization_cluster_counts.png', width=25, height=12, units='cm')

clusters.extended <- extend_clusters(clusters[cM_threshold==0.2])
fwrite(clusters.extended, file='data/clusters_extended.tsv', quote=F, row.names=F, col.names=T, sep='\t')

fwrite(clusters, file='data/clusters.tsv', quote=F, row.names=F, col.names=T, sep='\t')

cluster_sizes <- clusters[, list('N_clusters'=.N, 'N_singles'=sum(nSNPs_in_cluster==1)), by=list(cM_threshold, GWAS)]
cluster_sizes.long <- melt(cluster_sizes, measure.vars=c('N_clusters','N_singles'), value.name='count')
