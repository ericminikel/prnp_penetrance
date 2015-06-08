#!/bin/bash

# PRNP CDS:     20:4679867-4680628
# PRNP CDS +10: 20:4679857-4680638

# change to working directory and load paths to VCF files and metadata into current environment
. private-paths.bash

# extract PRNP from VCFs
tabix -h $exacvcf   20:4679857-4680638 > ./genotypes.vcf
tabix -h $exacsites 20:4679857-4680638 > ./annotations.vcf

# run Python script to grab everything of interest from these VCFs
bsub -q priority -W 4:00 -P $RANDOM -M 8000000 -J extract -o extract.o -e extract.e "./extract-variant-individuals.py -g genotypes.vcf.gz -a annotations.vcf -m $metadata --mingq 10 --minad 3 --minab .2 --maxmaf .001 > extract-variants.log 2> extract-variants.err"

# generate screenshots
mkdir -p igv
bsub -q priority -W 10:00 -P $RANDOM -M 4000000 -J igv -o igv/igvjob.out -e igv/igvjob.err \
    /home/unix/mlek/bin/IGV_plotter/IGV_plotter -LOCUS_SAMPLES ./sample_locus.txt -SAMPLES_TO_BAMS ./samples.txt -SNAPSHOT_DIR ./igv

# get 1kg population information
python src/get_populations.py --pcs data_nosync/pca_exac_all.csv --weights data_nosync/pc_weights.txt --ped_1kg_path data_nosync/integrated_call_samples.20130502.ALL.ped --samples data_nosync/list_60_5k_no_space.tsv > data_nosync/exac_60706_pops.tsv
# get codon 129 information
summarize_genotypes.py --vcf $exac63kgenos --samples indivs_with_path_alleles.txt --chrom 20 --pos 4680251 > path_allele_codon129.txt

# get AN for relevant populations
cat pops_60.5k.txt | awk -v FS="\t" '$2 == "JPT" {print $1}' > jpt.tsv
cat pops_60.5k.txt | awk -v FS="\t" '$2 == "TSI" {print $1}' > tsi.tsv
summarize_genotypes.py --vcf $exac_fullvcf --samples jpt.tsv --chrom 20 --pos 4680561 > jpt_m232r.tsv
summarize_genotypes.py --vcf $exac_fullvcf --samples tsi.tsv --chrom 20 --pos 4680494 > tsi_v210i.tsv
summarize_genotypes.py --vcf $exac_fullvcf --samples jpt.tsv --chrom 20 --pos 4680404 > jpt_v180i.tsv
cat jpt_m232r.tsv | grep -v "\\./\\." | wc -l # number of individuals with calls - 100% here
cat tsi_v210i.tsv | grep -v "\\./\\." | wc -l # number of individuals with calls 100%
cat jpt_v180i.tsv | grep -v "\\./\\." | wc -l # number of individuals with calls 100%

