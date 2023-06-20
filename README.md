# AGAAT
* **A**utomated **G**enotyping **A**rray **A**nalysis **T**ool (AGAAT) for Genotyping Array analysis.

## Dependencies

* Python version = 3.8.8

* OS = Ubuntu 20.04.4 

* Illumina Array Analysis Platform Genotyping Command Line Interface (iaap-cli) (v1.1)

* Affymetrix Analysis Power Tools (APT) (apt_2.11.6_linux_64_bit_x86_binaries.zip)

* PLINK (v1.90b6.21)

* samtools (1.12)

* bcftools (1.16)

* bcftools +gtc2vcf and +affy2vcf (version 2022-12-21 https://github.com/freeseek/gtc2vcf)

* gzip 

* pandas (python library, version=1.3.1)



### Required Files
* Human reference genome (hg19.fasta)

* bpm_manifest_file (--bpm GSA-24v3-0_A1.bpm, Infinium Global Screening Array Manifest File - BPM Format)

* csv_manifest_file (--csv GSA-24v3-0_A1.csv , Infinium Global Screening Array Manifest File - CSV Format)

* egt_cluster_file (--egt GSA-24v3-0_A1_ClusterFile.egt, Infinium Global Screening Array  Cluster File)

* idat_files (directory having all idat files)

* cel_files (directory having all cel files)

* Axiom APMRA Analysis and Annotation files

* phenotype(--pheno pheno1.txt) file with family and individual IDs of case samples in the first two columns

* dbsnp database vcf file (--dbsnp_common "common_all_20180418.vcf") with all the common variants

* gene list (--gene_list "gene_list.txt") for selecting markers belonging to the specified category of genes
  
* gene_bed file (bed file with genomic coordinates of genes), can be downloaded from UCSC table browser as UCSC_canonical.bed

```

  The file UCSC_canonical.bed looks like:

#hg19.knownCanonical.chrom      hg19.knownCanonical.chromStart  hg19.knownCanonical.chromEnd    hg19.knownCanonical.transcript  hg19.kgXref.geneSymbol
chr1    11873   14409   uc010nxq.1      DDX11L1
chr1    14361   19759   uc009viu.3      WASH7P
chr1    14406   29370   uc009viw.2      WASH7P
chr1    34610   36081   uc001aak.3      FAM138F
chr1    69090   70008   uc001aal.1      OR4F5
chr1    134772  140566  uc021oeg.2      LOC729737
chr1    321083  321115  uc001aaq.2      DQ597235
chr1    321145  321207  uc001aar.2      DQ599768
chr1    322036  326938  uc009vjk.2      LOC100133331

```

## Scripts

**agaat** has scripts for various genotyping array data analysis

The following scripts are available with this repo:
```
gsa_pipeline.sh

--Bash script pipeline for converting raw idat files into vcf files, followed by case-control association analysis using PLINK
```

```
merge_plink.sh

--Bash script pipeline for merging separate set of vcf files with an existing binary fileset of PLINK
```

```
apmra_pipe.sh

--Bash script pipeline for converting raw cel files of APMRA into vcf files, followed by case-control association analysis using PLINK
```

```
candidate_gene.py

--Python script for multiple testing correction using genetic variants belonging to a category of genes
```

```
ld_blocks.py

---Python script for multiple testing correction using haplotype blocks estimated by PLINK
```

```
candidate_ld.py

---Python script for candidate gene subset analysis and multiple testing correcting using haplotype blocks
```

```
common_var.py

---Python script for multiple testing correction using common variants belonging to a list of genes 
```

```
commonvar_ld.py

---Python script for candidate-gene, common-variant association analysis and multiple testing correcting using haplotype blocks
```

```
add_controls.py

---Python script for adding genotype counts from 1000 genome phase-3 populations to control data
```

## Usage

### Running the software

```
git clone https://github.com/anilprakash94/agaat.git agaat

cd agaat

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
	 -R,--ref_code <string>   Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
	 -p,--pheno <file>        Text file with family and individual IDs of case samples in the first two columns.
	 -t,--thresh <float>      p-value threshold of Hardy-Weinberg equilibrium test for filtering out variants.

```
--pheno file example: "pheno1.txt"

```
FID IID
0 205247220003_R01C01
0 205247220003_R01C02
0 205247220003_R02C01
0 205247220003_R02C02
```

```
bash merge_plink.sh -h

usage: bash merge_plink.sh [OPTIONS]
	 -h,--help                  Prints this message.
	 -s,--src2_vcf <directory>  Path of directory with vcf files to be added.
	 -R,--ref_code <string>     Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
	 -S,--src1 <file_prefix>    Prefix of existing PLINK binary fileset
	 -p,--pheno <file>          Phenotype text file of the new fileset with family and individual IDs of case samples in the first two columns.
         -t,--thresh <float>        p-value threshold of Hardy-Weinberg equilibrium test for filtering out variants.

```

```
bash apmra_pipe.sh -h

usage: bash apmra_pipe.sh [OPTIONS]
	 -h,--help                Prints this message.
	 -i,--cel_dir <dir>       Directory with input .cel files.
	 -a,--an_zip <file>       Compressed zip file having all the analysis files.
	 -x,--an_dir <dir>        Directory with all the analysis files.
	 -z,--annot_zip <file>    Compressed zip file having all the annotation files.
	 -d,--dqc_xml <file>      XML having paramaters for DQC value generation.
	 -c,--cr_xml <file>       XML having paramaters for QC call rates.
	 -s,--summ_xml <file>     XML having paramaters for summary intensity signals.
	 -n,--cnv_xml <file>      XML having paramaters for CNV analysis.
	 -g,--geno_xml <file>     XML having paramaters for genotype calls.
	 -o,--annot_csv <file>    Annotation file in CSV format.
	 -r,--ref_gen <file>      Human reference genome fasta file.
	 -v,--out_vcf <file>      Output VCF file generated from .txt files.
	 -p,--pheno <file>        Text file (space-delimited) with family and individual IDs of case samples in the first two columns.
	 -R,--ref_code <string>   Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
	 -t,--thresh <float>      P-value threshold of Hardy-Weinberg equilibrium test for filtering out variants.

```

```
python3 candidate_gene.py -h

usage: candidate_gene.py [-h] [--gene_list GENE_LIST] [--gene_bed GENE_BED]
                         [--plink_assoc PLINK_ASSOC] [--out_file OUT_FILE]

Multiple testing correction using genetic variants belonging to a category of
genes

optional arguments:
  -h, --help            show this help message and exit
  --gene_list GENE_LIST
                        input file having list of genes
  --gene_bed GENE_BED   bed file with genomic coordinates of genes
  --plink_assoc PLINK_ASSOC
                        plink association output file
  --out_file OUT_FILE   output file with gene-list based associations


```
--gene_bed can be extracted from UCSC table browser

--gene_list file example: "gene_list.txt"

```
AAK1
ABHD17A
ABHD17B
ABHD17C
ABHD6
ABI1
ABI2
ABI3
ABL1
ABL2
ABLIM3
```

Estimate haplotype blocks from the binary fileset using plink

```
usage : plink --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 500 --bfile bin_prefix

--bfile bin_prefix is the prefix of PLINK binary fileset

--blocks-max-kb is the maximum kilobase limit, variant pairs within this limit are only considered

```

The command above will generate "plink.blocks" file which can be used for LD-based multiple testing correction.

```

Multiple testing correction can be done using independent markers identified from haplotype blocks

python3 ld_blocks.py -h

usage: ld_blocks.py [-h] [--block_file BLOCK_FILE] [--assoc_adj ASSOC_ADJ]
                    [--out_file OUT_FILE]

Multiple testing correction using haplotype blocks

optional arguments:
  -h, --help            show this help message and exit
  --block_file BLOCK_FILE
                        plink --blocks output file with haplotype blocks and
                        variant IDs
  --assoc_adj ASSOC_ADJ
                        plink file with adjusted associations
  --out_file OUT_FILE   output file with ld-block based adjusted associations

```

```
Association analysis can be restricted to a candidate gene subset, followed by multiple testing correcting using haplotype blocks

python3 candidate_ld.py -h
usage: candidate_ld.py [-h] [--block_file BLOCK_FILE] [--gene_list GENE_LIST] [--gene_bed GENE_BED]
                       [--plink_assoc PLINK_ASSOC] [--out_file OUT_FILE]

Candidate gene association analysis and multiple testing correcting using haplotype blocks

optional arguments:
  -h, --help            show this help message and exit
  --block_file BLOCK_FILE
                        plink --blocks output file with haplotype blocks and variant IDs
  --gene_list GENE_LIST
                        input file having list of genes
  --gene_bed GENE_BED   bed file with genomic coordinates of genes
  --plink_assoc PLINK_ASSOC
                        plink association output file
  --out_file OUT_FILE   output file with gene-list based associations

  --out_file OUT_FILE   output file with ld-block based adjusted associations

```

```
python3 common_var.py -h
usage: common_var.py [-h] [--gene_list GENE_LIST] [--dbsnp_common DBSNP_COMMON] [--plink_adj PLINK_ADJ]
                     [--out_file OUT_FILE]

Multiple testing correction using common variants belonging to a list of genes

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

```
python3 commonvar_ld.py -h
usage: commonvar_ld.py [-h] [--block_file BLOCK_FILE] [--gene_list GENE_LIST] [--dbsnp_common DBSNP_COMMON]
                       [--plink_adj PLINK_ADJ] [--out_file OUT_FILE]

Candidate-gene common-variant association analysis and multiple testing correcting using haplotype blocks

optional arguments:
  -h, --help            show this help message and exit
  --block_file BLOCK_FILE
                        plink --blocks output file with haplotype blocks and variant IDs
  --gene_list GENE_LIST
                        input file having list of genes
  --dbsnp_common DBSNP_COMMON
                        dbsnp vcf file with common variants
  --plink_adj PLINK_ADJ
                        plink file with adjusted associations
  --out_file OUT_FILE   output file with ld-block based adjusted associations

```

```
	
Add genotype counts to control data followed by association test

Run in bash the following command for repeated runs until completion (The program may terminate due to connectivity issues, "until" commmand will resume the runs automatically)
#until python3 add_controls.py; do :; done

#plink --freqx --bfile final_merge --allow-no-sex --filter-controls
The above command produces "plink.frqx" genotype count report file for control samples

The input file, "plink.assoc", can be sorted based on raw p-value and a new text file with the top variants can be used as input to reduce time requirements

Output file from this code, 'assoc_1kgen_controls' can be used as an input instead of 'plink.assoc.adjusted' in other Python scripts with this repository
	
python3 add_controls.py -h
usage: add_controls.py [-h] [--assoc_num ASSOC_NUM] [--assoc_file ASSOC_FILE]
                       [--frqx_file FRQX_FILE] [--thresh THRESH]
                       [--case_num CASE_NUM] [--out_file OUT_FILE]
                       [--pop_codes POP_CODES] [--var_file VAR_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --assoc_num ASSOC_NUM
                        number of variants in the input file, give value as 0
                        if entire plink.assoc file is used as input
  --assoc_file ASSOC_FILE
                        plink association output file
  --frqx_file FRQX_FILE
                        plink genotype count report file
  --thresh THRESH       p-value threshold of Hardy-Weinberg equilibrium test
                        for filtering out variants
  --case_num CASE_NUM   number of case samples
  --out_file OUT_FILE   output file with hwe filtered and adjusted
                        associations
  --pop_codes POP_CODES
                        codes of 1000 genome phase-3 populations from which
                        control genotype counts will be added, comma-separated
  --var_file VAR_FILE   text file having variant data in which genotype counts
                        are already added


```
