#!/usr/bin/env bash


#!/usr/bin/env bash
INFILE=${1}
OUTFILE="${INFILE%.tsv.gz}.tab.bgz"

module load samtools


# HEADER="Chr\tBP\tSNP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tSE\tp\t"

zcat ${INFILE} | \
    awk 'BEGIN{OFS="\t"} NR>1 {print $2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' | \
    sort -nk2 | \
    bgzip -c > ${OUTFILE}

# build tabix index
tabix -f -b 2 -e 2 ${OUTFILE}

# tabix Hippocampus-EUR_18.tab.bgz 18:8691581-8691993
