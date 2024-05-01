#!/usr/bin/env bash
#SBATCH -J samtools
#SBATCH --partition=medium
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

SAM=$1

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate sam_bcf_tools_env

Short=$(basename $SAM .sam)

samtools view -@ 8 -bS $SAM -o $Short.bam
samtools sort -@ 8 $Short.bam -o "$Short"_sorted.bam
samtools index -@ 8 "$Short"_sorted.bam

# CHECK SORTED BAM CREATED
if [ -s "$Short"_sorted.bam ]; then
    rm $SAM
    rm $Short.bam
else
    echo "ERROR: "$Short"_sorted.bam is empty or does not exist."
fi

