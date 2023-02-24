#!/bin/bash

args=$(getopt --options s:R:S:p:h --longoptions src2_vcf:ref_code:src1:pheno:help -- "$@")

eval set -- "$args"


src2_vcf="/media/hmg/twocases_vcf"
ref_code="hg19"
src1="source1"
pheno="pheno2.txt" #phenotype file of the new dataset with family and individual IDs of case samples in the first two columns


while :
do
  case "$1" in
    -s | --src2_vcf )
      src2_vcf="$2"
      shift 2
      ;;
    -R | --ref_code )
      ref_code="$2"
      shift 2
      ;;
    -S | --src1 )
      src1="$2"
      shift 2
      ;;
    -p | --pheno )
      pheno="$2"
      shift 2
      ;;
    -h | --help)
      echo "For merging a separate set of vcf files with an existing binary fileset of PLINK.
      merge_plink Command-line Arguments
	  ================================
	  usage: bash merge_plink.sh [OPTIONS]
	 -h,--help                  Prints this message.
	 -s,--src2_vcf <directory>  Path of directory with vcf files to be added.
	 -R,--ref_code <string>     Human reference genome build code for PLINK : 'b36'/'hg18', 'b37'/'hg19', 'b38'/'hg38'.
	 -S,--src1 <file_prefix>    Prefix of existing PLINK binary fileset
	 -p,--pheno <file>          Phenotype text file of the new fileset with family and individual IDs of case samples in the first two columns."
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

if [[ "$src2_vcf" != */ ]]
then
  src2_vcf=$src2_vcf"/"
fi

#index compressed files
name=$(ls $src2_vcf | grep ".vcf.gz$")

for vcf_file in `echo "$name"`; do
bcftools index $src2_vcf$vcf_file
done

#create list of all vcf files of the second source
find $src2_vcf | grep -h ".vcf.gz$" > $src2_vcf"list.txt"

bcftools merge -l $src2_vcf"list.txt" -Oz -o src2_files.vcf.gz

gzip -d src2_files.vcf.gz


plink --vcf src2_files.vcf --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x $ref_code no-fail --allow-no-sex --make-bed --make-pheno $pheno '*' --out source2


plink --bfile $src1 --bmerge source2 --indiv-sort 0 -out merge

file="merge.missnp"

if test -f "$file"; then

    plink --bfile $src1 --keep-allele-order --const-fid --allow-extra-chr 0 --exclude merge.missnp --allow-no-sex --make-bed --out source1_tmp

    plink --bfile source2 --keep-allele-order --const-fid --allow-extra-chr 0 --exclude merge.missnp --allow-no-sex --make-bed --out source2_tmp

    plink --bfile source1_tmp --bmerge source2_tmp --indiv-sort 0 --allow-no-sex --make-bed --out final_merge

    rm source1_tmp.*

    rm source2_tmp.*

    plink --assoc counts --adjust --bfile final_merge --allow-no-sex --geno --mind
    
else
    
    plink --assoc counts --adjust --bfile merge --allow-no-sex --geno --mind
    
fi

