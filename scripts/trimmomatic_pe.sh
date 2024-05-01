#!/usr/bin/env bash
#SBATCH -J trimmomatic
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4

# F reads = $1 
# R reads = $2

ln -s $1
ln -s $2

file1=$(basename $1)
file2=$(basename $2)
fileshort=$(basename $1 _1.fastq.gz)

trimmomatic PE -threads 16 -phred33 $file1 $file2 \
    "$fileshort"_F_paired.fastq.gz "$fileshort"_F_unpaired.fastq.gz \
    "$fileshort"_R_paired.fastq.gz "$fileshort"_R_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 HEADCROP:10 MINLEN:80

rm $file1 $file2 "$fileshort"_F_unpaired.fastq.gz "$fileshort"_R_unpaired.fastq.gz

# CHECK OUTPUT CREATED
if [ -s "$fileshort"_F_paired.fastq.gz && -s "$fileshort"_R_paired.fastq.gz ]; then
    rm $file1 
    rm $file2 
    rm "$fileshort"_F_unpaired.fastq.gz 
    rm "$fileshort"_R_unpaired.fastq.gz
else
    echo "ERROR: Output is empty or does not exist."
fi
