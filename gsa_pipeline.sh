#!/bin/bash

args=$(getopt --options x:b:c:e:i:g:r:R:p:h --longoptions iaap:bpm:csv:egt:idat:gtc:ref:ref_code:pheno:help -- "$@")

eval set -- "$args"

iaap="/home/hmg/ngs/illumina_gsa_analysis/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/iaap-cli"
bpm="GSA-24v3-0_A1.bpm"
csv="GSA-24v3-0_A1.csv"
egt="GSA-24v3-0_A1_ClusterFile.egt"
idat="idat_files"
gtc="gtc_files"
ref="hg19.fasta"
ref_code="hg19" #human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'
pheno="pheno.txt" #pheno.txt file with family and individual IDs of case samples in the first two columns

while :
do
  case "$1" in
    -x | --iaap )
      iaap="$2"
      shift 2
      ;;
    -b | --bpm )
      bpm="$2"
      shift 2
      ;;
    -c | --csv )
      csv="$2"
      shift 2
      ;;
    -e | --egt )
      egt="$2"
      shift 2
      ;;
    -i | --idat )
      idat="$2"
      shift 2
      ;;
    -g | --gtc )
      gtc="$2"
      shift 2
      ;;
    -r | --ref )
      ref="$2"
      shift 2
      ;;
    -R | --ref_code )
      ref_code="$2"
      shift 2
      ;;
    -p | --pheno )
      pheno="$2"
      shift 2
      ;;
    -h | --help)
      echo "Converts illumina gsa raw idat files into vcf files, followed by case-control association analysis using PLINK.
      	gsa_pipeline Command-line Arguments
	================================
	usage: bash gsa_pipeline.sh [OPTIONS]
	 -h,--help                Prints this message.
	 -x,--iaap <executable>   Illumina Array Analysis Platform Genotyping Command Line Interface Executable.
	 -b,--bpm <file>          Infinium Global Screening Array Manifest File - BPM Format.
	 -c,--csv <file>          Infinium Global Screening Array Manifest File - CSV Format.
	 -e,--egt <file>          Infinium Global Screening Array  Cluster File.
	 -i,--idat <directory>    Directory with all idat files.
	 -g,--gtc <directory>     Directory for saving gtc output.
	 -r,--ref <file>          Human reference genome fasta file.
	 -R,--ref_code <string>  Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
	 -p,--pheno <file>        Text file with family and individual IDs of case samples in the first two columns."
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Invalid option: $1"
      ;;
  esac
done


echo "Converting idat files to gtc"

LR_ICU_VERSION_OVERRIDE="$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" $iaap gencall $bpm $egt $gtc --idat-folder $idat --output-gtc --gender-estimate-call-rate-threshold -0.1

echo "Indexing reference genome"
samtools faidx $ref

echo "Converting gtc files into a single vcf file"

bcftools +gtc2vcf \
  --no-version -Ov \
  --bpm $bpm \
  --csv $csv \
  --egt $egt \
  --gtcs $gtc \
  --fasta-ref $ref \
  --extra all_files.tsv | \
  bcftools sort -Ov -T ./bcftools. | \
  bcftools norm --no-version -Oz -c x -f $ref > all_files.vcf.gz


echo "Extracting all_files.vcf from all_files.vcf.gz"

gzip -d all_files.vcf.gz

echo "Creating plink binary files"
plink --vcf all_files.vcf --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x $ref_code no-fail --allow-no-sex --make-bed --make-pheno $pheno '*' --out MyVars

echo "Plink asscoiation test"

plink --assoc counts --adjust --bfile MyVars --allow-no-sex --geno --mind

