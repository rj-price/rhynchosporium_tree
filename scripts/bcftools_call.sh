#!/usr/bin/env bash
#SBATCH -J bcftools
#SBATCH --partition=long
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

Reference=$1
BamList=$2

Prefix=$(basename $BamList .list)

bcftools mpileup --threads 8 -Ou -m 10 -f $Reference -b $BamList | \
    bcftools call --threads 8 -cv --ploidy 1 -o $Prefix.vcf