#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=2

# reads = $1

mkdir fastqc
fastqc $1 -t 2 -o fastqc
