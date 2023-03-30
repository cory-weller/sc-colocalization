#!/usr/bin/env bash

TISSUE=${1}
SMR='/data/CARD/projects/singlecell_humanbrain/coloc/smr-1.3.1'
OUTDIR="/data/CARD/projects/singlecell_humanbrain/coloc_20230322/"
TMPDIR="/lscratch/${SLURM_JOB_ID}"
eQTL_DIR='/data/CARD/projects/omicSynth/SMR_omics/eQTLs'
FILESTEM="2020-05-26-${TISSUE}"
BESD="${eQTL_DIR}/${FILESTEM}/${FILESTEM}-${CHR}-SMR-besd"

# bash ~/sc-colocalization/submitSMR.sh Basalganglia-EUR
# bash ~/sc-colocalization/submitSMR.sh Cerebellum-EUR
# bash ~/sc-colocalization/submitSMR.sh Cortex-AFR
# bash ~/sc-colocalization/submitSMR.sh Cortex-EUR
# bash ~/sc-colocalization/submitSMR.sh Hippocampus-EUR
# bash ~/sc-colocalization/submitSMR.sh Spinalcord-EUR


run_smr() {
    local chr=${1}
    local besd="${eQTL_DIR}/${FILESTEM}/${FILESTEM}-${chr}-SMR-besd"
    echo "Running SMR for ${TISSUE} chr ${chr}"

    # print ALL eQTL results regardless of p-value
    ${SMR} \
        --descriptive-cis \
        --beqtl-summary \
        ${besd} \
        --query 1 \
        --out ${TISSUE}_${chr}
    
    echo "gzipping output to ${OUTDIR}/${TISSUE}_${chr}.tsv.gz"
    # cleanup
    gzip -c ${TISSUE}_${chr}.txt > ${OUTDIR}/${TISSUE}_${chr}.tsv.gz
    rm ${TISSUE}_${chr}.txt
}


# RUN
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR} && cd ${TMPDIR}

echo "Tissue type ${TISSUE}"

for CHR in $(seq 1 22); do
    run_smr ${CHR}
done

echo "all done"