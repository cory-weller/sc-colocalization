#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=4)
library(ggplot2)

# read in formatted summary statistics - gr38

# AD_fn <- 'data/AD_Bellenguez.formatted.tsv.gz'
# LBD_fn <- 'data/LBD_Chia.formatted.tsv.gz'
# PD_fn <- 'data/PD_Nalls.formatted.tsv.gz'

# basalganglia_fn <- 'data/Basalganglia-EUR-out.cis.summary.txt'
# cerebellum_fn <- 'data/Cerebellum-out.cis.summary.txt'
# cortex_fn <- 'data/Cortex-EUR-out.cis.summary.txt'
# hippocampus_fn <- 'data/Hippocampus-out.cis.summary.txt'
# spinalcord_fn <- 'data/Spinalcord-out.cis.summary.txt'

diagnoses <- list(c('AD','data/AD_Bellenguez.formatted.tsv.gz'),
                    c('LBD','data/LBD_Chia.formatted.tsv.gz'),
                    c('PD','data/PD_Nalls.formatted.tsv.gz')
)

celltypes <- list(c('basal_ganglia','data/Basalganglia-EUR-out.cis.summary.txt'),
                c('cerebellum','data/Cerebellum-out.cis.summary.txt'),
                c('cortex','data/Cortex-EUR-out.cis.summary.txt'),
                c('hippocampus','data/Hippocampus-out.cis.summary.txt'),
                c('spinal_cord','data/Spinalcord-out.cis.summary.txt')
)

get_diagnosis <- function(diagnosis_fn, diagnosis_name) {
    dat <- fread(diagnosis_fn)
    dat[, c('SNP','A1','A2') := NULL]
    setnames(dat, c('BETA','FRQ','SE','P'), c('diagnosis_b','diagnosis_Freq', 'diagnosis_se', 'diagnosis_p'))
    setcolorder(dat, c('CHR','BP','diagnosis_Freq','diagnosis_b','diagnosis_se','diagnosis_p'))
    dat[, 'diagnosis' := diagnosis_name]
    return(dat[])
}

get_eQTL <- function(celltype_fn, celltype_name) {
    dat <- fread(celltype_fn)
    dat <- dat[SNP_Chr != 'SNP_Chr']
    setnames(dat, 'SNP_Chr', 'CHR')
    dat[, CHR := as.numeric(CHR)]
    dat[, BP := as.numeric(BP)]
    dat[, p := as.numeric(p)]
    dat[, se := as.numeric(se)]
    dat[, b := as.numeric(b)]
    dat[, Freq := as.numeric(Freq)]
    dat[,c('TopSNP','A1','A2','cis_Chr','cis_left','cis_right','Orientation','Probe_Chr','Probe', 'Probe_bp', 'Gene') := NULL]
    setnames(dat, c('Freq','b','se','p'), c('eQTL_Freq','eQTL_b','eQTL_se','eQTL_p'))
    setkey(dat, CHR, BP)
    dat[, 'celltype' := celltype_name]
    return(dat[])
}


dat <- foreach(i=diagnoses, .combine='rbind') %do% {
    diagnosis_type <- i[1]
    diagnosis_filename <- i[2]
    dt.diagnosis <- get_diagnosis(diagnosis_filename, diagnosis_type)
    foreach(j=celltypes, .combine='rbind') %do% {
        celltype_name <- j[1]
        celltype_filename <- j[2]
        dt.eQTL <- get_eQTL(celltype_filename, celltype_name)
        dt.merge <- merge(dt.diagnosis, dt.eQTL)
        return(dt.merge[])
    }
}



dat[, eQTL_hit := ifelse(eQTL_p < 0.05, TRUE, FALSE)]
dat[, diagnosis_hit := ifelse(diagnosis_p < 0.05, TRUE, FALSE)]



chrLengths <- data.table(CHR=1:22, bp=
c('248956422','242193529','198295559','190214555','181538259','170805979','159345973',
'145138636','138394717','133797422','135086622','133275309','114364328','107043718',
'101991189','90338345','83257441','80373285','58617616','64444167','46709983','50818468'))

chrLengths[, cumulative := shift(cumsum(bp), type='lag', fill=0)]
chrLengths[, bp := NULL]

dat  <- merge(dat, chrLengths, by.x='CHR', by.y='CHR')
fwrite(dat, file='merged_hits.tsv', quote=F, row.names=F, col.names=T, sep='\t')

dat <- dat[diagnosis_hit==TRUE & eQTL_hit==TRUE][order(celltype, diagnosis)]

dat[, grp := paste0(celltype, '_', diagnosis)]



dat[, cumulativeBP := BP + cumulative]
windowWidth <- 5e6
dat[, window := floor(cumulativeBP / windowWidth)]
dat.ag <- dat[, .N, by=list(grp, window)]

allcombos <- CJ('grp'=unique(dat.ag$grp), 'window'=0:max(dat.ag$window))
dat.ag.all <- merge(dat.ag, allcombos, all=TRUE)
dat.ag.all[is.na(N), N:=0]

g <- ggplot(data=dat.ag.all, aes(x=window, y=grp, fill=N)) +
geom_tile() +
scale_fill_viridis() +
labs(x='500 kb window', y='') +
theme_few() 

ggsave(g, file='heatmap.png', width=60, height=8, units='cm')