#!/usr/bin/env bash
#SBATCH -J raxml
#SBATCH --partition=long
#SBATCH --mem=30G
#SBATCH --cpus-per-task=12

PhylipFile=$1

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate sam_bcf_tools_env

# For just one ML tree
raxml-ng --msa $PhylipFile --model GTR+ASC_LEWIS 

#raxml-ng --threads 20 --workers 10 --all --msa $PhylipFile --model GTR+ASC_LEWIS --tree pars{10} --bs-trees 100

# --all means do a full analysis: estimate best ML tree, then do bootstraps, then use bootstraps to plot support on best tree
# --msa means multiple-sequence alignment
# --model indicates evolutionary model, here we use General Time Reversible plus the Lewis ascertainment bias correction (takes into account the fact that we are only using variable sites)
# --bs-trees indicates the number of bootstraps you'd like to do