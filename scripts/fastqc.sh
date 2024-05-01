#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=2

Reads=$1

mkdir fastqc
fastqc $Reads -t 2 -o fastqc
