#!/usr/bin/env bash
#SBATCH -J bwa-mem
#SBATCH --partition=medium
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

Genome=$1
F_Reads=$2
R_Reads=$3

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate gatk

Prefix=$(basename $2 _F.fastq.gz)

bwa mem -t 16 $Genome $F_Reads $R_Reads > $Prefix.sam
