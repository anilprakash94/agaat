#!/bin/bash

#SNP-specific priors are recommended for samples <= 96

args=$(getopt --options i:a:x:z:d:c:s:n:g:o:r:v:p:R:t:h --longoptions cel_dir:,an_zip:,an_dir:,annot_zip:,dqc_xml:,cr_xml:,summ_xml:,cnv_xml:,geno_xml:,annot_csv:,ref_gen:,out_vcf:,pheno:,ref_code:,thresh:,help -- "$@")

eval set -- "$args"


cel_dir="/home/affymetrix/test/cel_files"
an_zip="/home/affymetrix/support_files/APMRA/TFS-Assets_LSG_Support-Files_Axiom_APMRA_Analysis.r3.zip"
an_dir="/home/affymetrix/analysis_files"
annot_zip="/home/affymetrix/support_files/APMRA/TFS-Assets_LSG_Support-Files_Axiom_APMRA.na35.r3.a2.annot.csv.zip"
dqc_xml="Axiom_APMRA.r3.apt-geno-qc.AxiomQC1.xml"
cr_xml="Axiom_APMRA_SNPSpecificPriors_Step1.r3.apt-genotype-axiom.AxiomGT1.apt2.xml"
summ_xml="Axiom_APMRA.r3.apt-genotype-axiom.AxiomCN_GT1.apt2.xml"
cnv_xml="Axiom_APMRA.r3.apt-copynumber-axiom-hmm.AxiomHMM.apt2.xml"
geno_xml="Axiom_APMRA_SNPSpecificPriors_Step2.r3.apt-genotype-axiom.AxiomGT1.apt2.xml"
annot_csv="/home/affymetrix/analysis_files/Axiom_APMRA.na35.r3.a2.annot.csv"
ref_gen="/home/hg19.fasta"
out_vcf="/home/affymetrix/test/cel_files/affy_files.vcf"
pheno="/home/affymetrix/test/cel_files/pheno.txt" 
ref_code="hg19"
thresh="1e-5"


while :
do
  case "$1" in
    -i | --cel_dir )
      cel_dir="${2%/}"
      shift 2
      ;;
    -a | --an_zip )
      an_zip="$2"
      shift 2
      ;;
    -x | --an_dir )
      an_dir="${2%/}"
      shift 2
      ;;
    -z | --annot_zip )
      annot_zip="$2"
      shift 2
      ;;
    -d | --dqc_xml )
      dqc_xml="$2"
      shift 2
      ;;
    -c | --cr_xml )
      cr_xml="$2"
      shift 2
      ;;
    -s | --summ_xml )
      summ_xml="$2"
      shift 2
      ;;
    -n | --cnv_xml )
      cnv_xml="$2"
      shift 2
      ;;
    -g | --geno_xml )
      geno_xml="$2"
      shift 2
      ;;
    -o | --annot_csv )
      annot_csv="$2"
      shift 2
      ;;
    -r | --ref_gen )
      ref_gen="$2"
      shift 2
      ;;
    -v | --out_vcf )
      out_vcf="$2"
      shift 2
      ;;
    -p | --pheno )
      pheno="$2"
      shift 2
      ;;
    -R | --ref_code )
      ref_code="$2"
      shift 2
      ;;
    -t | --thresh )
      thresh="$2"
      shift 2
      ;;
    -h | --help)
      echo "Converts Axiom PMRA array raw .cel files into vcf files, followed by case-control association analysis using PLINK.
      	apmra_pipe Command-line Arguments
	================================
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
	 -t,--thresh <float>      P-value threshold of Hardy-Weinberg equilibrium test for filtering out variants."
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




echo -e "\nWriting .cel file names into text file with 'cel_files' as header\n"

(echo cel_files; \ls -1 "$cel_dir/"*.CEL ) > $cel_dir"/cel_list1.txt"


echo -e "\nExtracting analysis files into folder\n"

unzip $an_zip -d $an_dir

echo -e "\nExtracting annotation files into folder\n"

unzip $annot_zip -d $an_dir


echo -e "\nGenerating sample 'DQC' values\n"

apt-geno-qc \
    --analysis-files-path $an_dir \
    --xml-file $an_dir"/"$dqc_xml \
    --cel-files $cel_dir"/cel_list1.txt" \
    --out-file $cel_dir"/qc.txt" \
    --log-file $cel_dir"/apt-geno-qc.log"


echo -e "\nFiltering .cel files with 'axiom_dishqc_DQC' value < 0.82 in the qc.txt file\n"


echo cel_files > $cel_dir"/cel_list2.txt"
while read p || [ -n "$p" ]; do
    file_name=`echo $p | awk '{print $1}'`
    dqc_value=`echo $p | awk '{print $18}'`
    if [[ `echo "$dqc_vale 0.82" | awk '{print ($1 > $2)}'` == 1 ]] ; then
        echo $cel_dir"/"$file_name >> $cel_dir"/cel_list2.txt"
    fi
done< <(grep -v '^#' $cel_dir"/qc.txt" | tail -n +2)



echo -e "\nGenerating sample qc call rates\n"

apt-genotype-axiom \
    --log-file $cel_dir"/apt-genotype-axiom.log" \
    --arg-file $an_dir"/"$cr_xml \
    --analysis-files-path $an_dir \
    --out-dir $cel_dir"/step1" \
    --dual-channel-normalization true \
    --table-output false \
    --cel-files $cel_dir"/cel_list2.txt"


echo -e "\nFiltering .cel files with QC call rate value < 97%\n"

echo cel_files > $cel_dir"/cel_list3.txt"
while read q || [ -n "$q" ]; do
    file_name=`echo $q | awk '{print $1}'`
    call_rate=`echo $q | awk '{print $3}'`
    if [[ `echo "$call_rate 97" | awk '{print ($1 > $2)}'` == 1 ]] ; then
        echo $cel_dir"/"$file_name >> $cel_dir"/cel_list3.txt"
    fi  
done< <(grep -v '^#' $cel_dir"/step1/AxiomGT1.report.txt" | tail -n +2)


#samples passing DQC and QC call rate filtering are only processed for genotyping

#For arrays with copy number-aware genotyping (CNAG) enabled, the following two extra steps should be executed. CNAG is enabled in APMRA. 

echo -e "\nGenerating summary intensity signals for all probesets\n"

apt-genotype-axiom \
    --analysis-files-path $an_dir \
    --arg-file $an_dir"/"$summ_xml \
    --cel-files $cel_dir"/cel_list3.txt" \
    --out-dir $cel_dir"/summary" \
    --log-file $cel_dir"/summary/apt2-axiom.log"


echo -e "\nCopy number analysis\n"

apt-copynumber-axiom-hmm  \
    --analysis-files-path $an_dir \
    --arg-file $an_dir"/"$cnv_xml \
    --summary-file $cel_dir"/summary/AxiomGT1.summary.a5" \
    --report-file $cel_dir"/summary/AxiomGT1.report.txt" \
    --out-dir $cel_dir"/cn" \
    --log-file $cel_dir"/cn/apt-copynumber-axiom.log"


echo -e "\nProducing genotype calls\n"

apt-genotype-axiom \
    --log-file $cel_dir"/step2/apt-genotype-axiom.log" \
    --arg-file $an_dir"/"$geno_xml \
    --analysis-files-path $an_dir \
    --out-dir $cel_dir"/step2" \
    --dual-channel-normalization true \
    --allele-summaries true \
    --genotyping-node:snp-posteriors-output true \
    --batch-folder $cel_dir"/step2" \
    --cel-files $cel_dir"/cel_list3.txt"



echo "Indexing reference genome"
samtools faidx $ref_gen


echo -e "\nConverting .txt output files from APT to VCF\n"
 
bcftools +affy2vcf \
    --csv $annot_csv \
    --fasta-ref $ref_gen \
    --calls $cel_dir"/step2/AxiomGT1.calls.txt" \
    --confidences $cel_dir"/step2/AxiomGT1.confidences.txt" \
    --summary $cel_dir"/step2/AxiomGT1.summary.txt" \
    --snp $cel_dir"/step2/AxiomGT1.snp-posteriors.txt" \
    --output $out_vcf


echo -e "\nFiltering qc fail samples from pheno file\n"

echo "FID IID" > $cel_dir"/pheno1.txt"
while read u || [ -n "$u" ]; do
    file_name=`echo $u | awk '{print $2}'`
    if grep -Rq "$file_name" $cel_dir"/cel_list3.txt" ; then
        echo $u >> $cel_dir"/pheno1.txt"
    fi  
done< <(tail -n +2 $pheno)


echo -e "\nCreating plink binary files\n"
plink --vcf $out_vcf --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x $ref_code no-fail --allow-no-sex --make-bed --make-pheno $cel_dir"/pheno1.txt" '*' --out $cel_dir"/source1"


echo -e "\nPlink asscoiation test\n"

plink --assoc counts --adjust --hwe $thresh 'midp' --bfile $cel_dir"/source1" --allow-no-sex --geno --mind



