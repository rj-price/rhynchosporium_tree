#!/usr/bin/env bash
#SBATCH -J bcftools
#SBATCH --partition=long
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

Reference=$1
BamList=$2

Prefix=$(basename $BamList .list)

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate sam_bcf_tools_env

bcftools mpileup --threads 8 -Ou -m 5 -f $Reference -b $BamList | \
    bcftools call --threads 8 -cv --ploidy 1 -o $Prefix.bcf