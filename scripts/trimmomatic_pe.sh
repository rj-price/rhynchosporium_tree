#!/usr/bin/env bash
#SBATCH -J trimmomatic
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4

# F reads = $1 
# R reads = $2

ln -s $1
ln -s $2

FWD=$(basename $1)
REV=$(basename $2)
Prefix=$(basename $1 _1.fastq.gz)

trimmomatic PE -threads 8 -phred33 $FWD $REV \
    "$Prefix"_F_paired.fastq.gz "$Prefix"_F_unpaired.fastq.gz \
    "$Prefix"_R_paired.fastq.gz "$Prefix"_R_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 HEADCROP:10 MINLEN:80

rm $FWD $REV "$Prefix"_F_unpaired.fastq.gz "$Prefix"_R_unpaired.fastq.gz

# CHECK OUTPUT CREATED
if [ -s "$Prefix"_F_paired.fastq.gz && -s "$Prefix"_R_paired.fastq.gz ]; then
    rm $FWD 
    rm $REV 
    rm "$Prefix"_F_unpaired.fastq.gz 
    rm "$Prefix"_R_unpaired.fastq.gz
else
    echo "ERROR: Output is empty or does not exist."
fi
