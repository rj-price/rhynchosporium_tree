# *Rhynchosporium* Phylogeny
Phylogenetic tree reconstruction using available *Rhynchosporium commune* data.

## Prepare data

### Download WGS data and genome from NCBI
Search for *Rhynchosporium* DNA data on SRA and download metadata (SraRunTable.txt). Using ENA, generate a script to download all accessions from BioProjects PRJNA327656<sup>[1](#ref1)</sup> and PRJNA419548, and download to cluster. 
```
bash scripts/ena-file-download-20230803-1315.sh
```

Download *R. commune* data and genome, and uncompress.
```
conda activate ncbi_datasets
datasets download genome accession GCA_900074885.1 --filename GCA_900074885.1.zip
```

### QC reads
QC raw sequencing data.
```
for file in *gz; do 
    sbatch ../scripts/fastqc.sh $file
done

multiqc fastqc/
```

Trim adapters and low quality sequence, and QC trimmed reads.
```
mkdir trimmed && cd trimmed

for file in ../reads/*_1.fastq.gz; do 
    file2=$(ls $file | sed s/"_1.fastq.gz"//g)
    sbatch ../scripts/trimmomatic_pe.sh $file "$file2"_2.fastq.gz
done

for file in *gz; do 
    sbatch ../scripts/fastqc.sh $file
done

multiqc fastqc/
```

Rename reads to "Run_geo_loc_name_country_Collection_Date_readDirection.fastq.gz"
```
python scripts/rename_reads.py SraRunTable.txt trimmed/ renamed_reads/
```

## References
<a id="ref1">1.</a> Mohd-Assaad N, McDonald BA, Croll D. Genome-Wide Detection of Genes Under Positive Selection in Worldwide Populations of the Barley Scald Pathogen. Genome Biol Evol. 2018 Apr 1;10(5):1315-1332. doi: 10.1093/gbe/evy087. PMID: 29722810; PMCID: [PMC5972619](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5972619).