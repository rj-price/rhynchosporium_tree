#!/usr/bin/env bash
#SBATCH -J bwa-mem
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

Genome=$1
F_Reads=$2
R_Reads=$3

Prefix=$(basename $2 _F.fastq.gz)

bwa mem -t 8 $Genome $F_Reads $R_Reads > $Prefix.sam
