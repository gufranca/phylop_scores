#!/bin/bash


# run Random Flanking simulation for intergenic (N = 100)
while read line; do
    echo $line | sed 's/ /\t/g' > tmp1
    python simulation_features.py -i tmp1 -b \
    ensembl71_protein_coding_exons.bed  -d bed_files/ \
    -rf -n 100 >> mirnas_7_12_inter_phylop_flank_rf100.txt
    rm -rf tmp1
done < mirnas_7_12_inter.bed


# run Random Intragenic simulation for intragenic miRNAs (N = 100)
while read line; do
    echo $line | sed 's/ /\t/g' > tmp1
    python simulation_features.py -i tmp1 -b \
    ensembl71_protein_coding_introns.bed  -d bed_files/ \
    -r -n 100 >> mirnas_7_12_intra_phylop_random_r100.txt
    rm -rf tmp1
done < mirnas_7_12_intra_less_exonic.bed


