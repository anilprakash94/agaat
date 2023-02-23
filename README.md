# gsa_pipeline
* Pipeline for illumina gsa analysis.

## Deependencies

* Python version = 3.8.8

* OS = Ubuntu 20.04.4 

* Illumina Array Analysis Platform Genotyping Command Line Interface (iaap-cli) (v1.1)

* PLINK (v1.90b6.21)

* samtools (1.12)

* bcftools (1.16)

* bcftools +gtc2vcf (version 2022-12-21 https://github.com/freeseek/gtc2vcf)

* gzip 

* pandas (python library, version=1.3.1)



### Required Files
* Human reference genome (hg19.fasta)

* bpm_manifest_file (Infinium Global Screening Array Manifest File - BPM Format)

* csv_manifest_file (Infinium Global Screening Array Manifest File - CSV Format)

* egt_cluster_file (Infinium Global Screening Array  Cluster File)

* idat_files (All raw idat files in a directory)

* pheno.txt file with family and individual IDs of case samples in the first two columns

## Scripts

**gsa pipeline** has scripts for various gsa data analysis

The following scripts are available with this repo:
```
gsa_pipeline.sh

--Bash script pipeline for converting raw idat files into vcf files
```
```
plink_categ_assoc.py

--Python script for multiple testing correction using genetic variants belonging to a category of genes
```


## Usage

### Running the model

```
git clone https://github.com/anilprakash94/gsa_pipeline.git gsa_pipeline

cd gsa_pipeline

```
Then, run the programs according to the requirements and instructions listed in README.md.

For example:


```
bash gsa_pipeline.sh -h

usage: bash gsa_pipeline.sh [OPTIONS]
     -h,--help                Prints this message.
     -x,--iaap <executable>   Illumina Array Analysis Platform Genotyping Command Line Interface Executable.
     -b,--bpm <file>          Infinium Global Screening Array Manifest File - BPM Format.
     -c,--csv <file>          Infinium Global Screening Array Manifest File - CSV Format.
     -e,--egt <file>          Infinium Global Screening Array  Cluster File.
     -i,--idat <directory>    Directory with all idat files.
     -g,--gtc <directory>     Directory for saving gtc output.
     -r,--ref <file>          Human reference genome fasta file.
     -rc,--ref_code <string>  Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
     -p,--pheno <file>        Text file with family and individual IDs of case samples in the first two columns.

```
```
python3 plink_categ_assoc.py -h

usage: plink_categ_assoc.py [-h] [--gene_list GENE_LIST]
                            [--dbsnp_common DBSNP_COMMON]
                            [--plink_adj PLINK_ADJ] [--out_file OUT_FILE]

Multiple testing correction using genetic variants belonging to a category of
genes

optional arguments:
  -h, --help            show this help message and exit
  --gene_list GENE_LIST
                        input file having list of genes
  --dbsnp_common DBSNP_COMMON
                        dbsnp vcf file with common variants
  --plink_adj PLINK_ADJ
                        plink file with adjusted associations
  --out_file OUT_FILE   output file with gene-list based adjusted associations

```
