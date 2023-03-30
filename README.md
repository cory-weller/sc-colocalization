# sc-colocalization

Generate full tissue eQTL GWAS files using [`submitSMR.sh`](src/submitSMR.sh). Importantly, the
script calls a `smr-1.3.1` binary with `--descriptive-cis` `--beqtl-summary` to generate output,
and `--query 1` to include
    
## Collecting significant tissue (eQTL GWAS) data
This loop iterates through (tissue) x (chromosomes 1:22) to query compressed `.besd` files. Records
passing the defined p-value significance threshold (in this case `1e-4`) are concatenated into a
single `${tissue}.signif.tsv.gz` file per tissue type, written to `data/TISSUE_eQTL`.

```bash
TISSUES=(Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EUR Hippocampus-EUR Spinalcord-EUR)
THRESHOLD='1e-4'
for tissue in ${TISSUES[@]}; do
    bash src/retrieve-smr.sh ${tissue} ${THRESHOLD}
done
```


## Collecting significant NDD GWAS data
This loop iterates through `${DATADIR}`, which currently contains multiple NDD GWAS files with
varied column order and ID. The function `filterPvalByR` envokes an `R` heredoc to read in the
complete GWAS, then output a consistent set of columns, with rows filtered by p-value significance.

```bash
filterPvalByR() {
Rscript - ${1} ${2} ${3} <<EOF
    #!/usr/bin/env Rscript
    library(data.table)
    args <- commandArgs(trailingOnly=TRUE)
    infile <- args[1]
    threshold <- as.numeric(args[2])
    outfile <- args[3]

    dat <- fread(infile)
    desired_cols <- c('SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P')
    dat <- dat[P < threshold, .SD, .SDcols=desired_cols]
    fwrite(dat, file=outfile, quote=F, row.names=F, col.names=T, sep='\t')
EOF
}

DATADIR='/gpfs/gsfs8/users/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats'
OUTDIR='/gpfs/gsfs9/users/wellerca/sc-colocalization/data/NDD'
THRESHOLD='5E-8'

for infile in ${DATADIR}/*; do
    outfile=$(basename ${infile%.formatted.tsv}.signif.tsv)
    filterPvalByR ${infile} ${THRESHOLD} ${OUTDIR}/${outfile}
done
```

## Defining clusters for colocalization analysis

## Running iterative colocalization


```bash
Rscript find_clusters.R
Rscript subset_eQTL_gwas.R
Rscript subset_diagnosis_gwas.R
Rscript run_coloc.R
```

Generating filtered 
```bash
filterPvalByAwk() {
    local infile=${1}
    local p_val_column=${2}
    local p_val_threshold=${3}
    zcat ${infile} | \
    awk -v col="${p_val_column}" \
        -v threshold="${p_val_threshold}" \
        'N == 1 || ($col + 0) < threshold'
}

datadir='/gpfs/gsfs9/users/CARD/projects/singlecell_humanbrain/coloc_20230322'


# Filtering eQTL GWAS files to only include SNPs with P < 1E-4
cd ~/sc-colocalization/data
TISSUES=(Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EUR Hippocampus-EUR Spinalcord-EUR)

for tissue in ${TISSUES[@]}; do
    for chr in $(seq 1 21); do
        echo "$tissue chr${chr}"
        filename=${datadir}/${tissue}_${chr}.tsv.gz
        filterPvalByAwk ${filename} 14 '1E-4' >> ${tissue}.allChr.signif.tsv
    done
done

gzip *.allChr.signif.tsv
eQTL_filename <- paste0('', TISSUE, '.allChr.signif.tsv.gz')
NDD_filename <- 
```

```bash


```



```

```




```
awk -v col='14' \
    -v threshold='1E-4' \
    '($col + 0) < threshold' \
    10_eQTL_all_SNPs.tsv | head

./smr-1.3.1 --descriptive-cis --beqtl-summary /data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Basalganglia-EUR/2020-05-26-Basalganglia-EUR-10-SMR-besd --out Basalganglia-EUR-10-
```
AD_Bellenguez.formatted.tsv.gz
PD_Nalls.formatted.tsv.gz
LBD_Chia.formatted.tsv.gz

# Get linkage map

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz


dirs:
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Basalganglia-EUR
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Cerebellum-EUR
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Cortex-AFR
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Cortex-EUR
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Hippocampus-EUR
/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Spinalcord-EUR

celltypes=(Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EUR Hippocampus-EUR Spinalcord-EUR)
diagnoses=(AD PD LBD)

For each celltype
    For each chromosome
        Submit job, request lscratch, cd to TMP
        generate full SNP table for celltype eQTL
        Start R
        load eQTL into R
        For each diagnosis (AD PD LBD)
            load diagnosis GWAS into R
            run coloc
            output results to permanent dir


GWAS files stored in
```bash
'/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/outfiles'
```

filted to significant GWAS hits with:

```R
#!/usr/bin/env Rscript

library(data.table)

disease_gwas_files <- list.files(path='/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats', full.names=T)
signif_threshold <- 5e-8

for(fn in disease_gwas_files) {
    dat <- fread(fn)
    base_filename <- basename(fn)
    gwas_name <- strsplit(base_filename, split='\\.formatted')[[1]][1]
    out_filename <- paste0(gwas_name, '.signif.tsv')

    dat <- dat[P <= signif_threshold]
    fwrite(dat, file=paste0('data/', out_filename), quote=F, row.names=F, col.names=T, sep='\t')
}
```

In R, load in bed file (cM) and significant GWAS
[, cM := to_cM(BP)]
[, nextSNP := shift(BP, size=1L)]
[, nextSNP_cM := to_cM(nextSNP)]
[, cM_between_SNPs := nextSNP_cM - cM]
[, tf := ifelse(cM_between_SNPs < 0.1, TRUE, FALSE)]    # TRUE if next SNP is within 0.1 cM, otherwise FALSE
[, grp := rleid(CHR,tf)]                                # Assigns group ID based on continuous runs of CHR and tf