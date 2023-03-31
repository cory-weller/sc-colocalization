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
```bash
RECOMBINATION_BED='data/genetic_map_hg38_withX.txt.gz'
wget -O ${RECOMBINATION_BED} https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
NDD_FILES=($(ls data/NDD/*signif.tsv))

for ndd_gwas in ${NDD_FILES[@]}; do
    Rscript src/find_clusters.R ${ndd_gwas} ${RECOMBINATION_BED} ${cM_THRESHOLD}
done
```

## Running iterative colocalization

```bash
eQTL_FILES=($(ls data/eQTL/*signif.tsv.gz))
NDD_FILES=($(ls data/NDD/*signif.tsv))

for eQTL_file in ${eQTL_FILES[@]}; do 
    for ndd_gwas in ${NDD_FILES[@]}; do
        clusterfile="data/NDD/clusters/$(basename ${ndd_gwas[0]%.signif.tsv}.clusters_chosen.tsv)"
        Rscript src/run-coloc.R ${eQTL_file} ${ndd_gwas} ${clusterfile}
    done
done
```

## Combining results files
```bash
OUTFILE='data/coloc/ALL_COLOC.tsv'
INFILES=($(ls data/coloc/*.tsv))

if [ -f "${OUTFILE}" ]; then rm ${OUTFILE}; fi

head -n 1 ${INFILES[0]} > ${OUTFILE}

for file in ${INFILES[@]}; do
    awk 'NR>1' ${file} >> ${OUTFILE}
done
```

```R
library(data.table)
dat <- fread('data/coloc/ALL_COLOC.tsv')
gids <- fread('data/coloc/GENEIDS.csv', col.names=c('PROBE','GENEID'))
dat <- merge(dat, unique(gids), by='PROBE')
setnames(dat, 'B', 'eQTL_B')

desired_cols <- c('CHR','BP','NDD','TISSUE','GENEID','SNP.PP.H4','eQTL_B','PROBE','PROBE_BP','SNP','A1','A2','FREQ')
dat.out <- dat[, .SD, .SDcols=desired_cols][order(CHR,BP,NDD,TISSUE)]
fwrite(dat.out, file='coloc_formatted.tsv', quote=F, row.names=F, col.names=T, sep='\t')
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


https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

