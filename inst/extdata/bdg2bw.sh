#!/bin/bash


F=$1
G=$2
bt=$3
bc=$4
bb=$5

cat ${F} | sed 's/ //g' | awk '{OFS="\t"; print $1,$2,$3,$4}' | ${bt} slop -i stdin -g ${G} -b 0 | ${bc} stdin ${G} ${F}.clip

LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip

${bb} ${F}.sort.clip ${G} ${F/bedGraph/bw}

mv ${F}.clip ${F}

rm -f ${F}.sort.clip ${F}.clip
