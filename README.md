# sc-colocalization

Generating full SNP results files from besd
```bash
module load SMR
smr --beqtl-summary $SMR_DIR/TEST_DATA/westra_eqtl_hg18 --query 5.0e-8 --snp rs123


https://hpc.nih.gov/apps/SMR.html#:~:text=SMR%20(Summary%2Dbased%20Mendelian%20Randomization,complex%20trait%20because%20of%20pleiotropy.

https://cran.r-project.org/web/packages/coloc/vignettes/a02_data.html

/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats

datadir='/data/CARD/projects/omicSynth/SMR_omics/eQTLs/2020-05-26-Cerebellum-EUR/'

smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-10-SMR-besd --query 1
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-11-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-12-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-13-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-14-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-15-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-16-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-17-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-18-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-19-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-1-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-20-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-21-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-22-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-2-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-3-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-4-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-5-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-6-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-7-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-8-SMR-besd
smr --descriptive-cis --beqtl-summary ${datadir}/2020-05-26-Cerebellum-EUR-9-SMR-besd
```

Subset eQTL GWAS tables to only those needed for colocalization, generating `data/eQTL_colocalization_SNPs.tsv` 
and `data/eQTL_colocalization_SNPs.tsv`
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

datadir='/gpfs/gsfs8/users/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats'
outdir='/gpfs/gsfs9/users/wellerca/sc-colocalization/data/NDD'
threshold='5E-8'

for infile in ${datadir}/*; do
    outfile=$(basename ${infile%.formatted.tsv}.signif.tsv)
    filterPvalByR ${infile} ${threshold} ${outdir}/${outfile}
done

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