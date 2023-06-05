import argparse, requests, re, os
import pandas as pd
import numpy as np
from scipy.stats import chisquare
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings

disable_warnings(InsecureRequestWarning)

#Run in bash the following command for repeated runs until completion
#until python3 add_controls.py; do :; done

#plink --freqx --bfile final_merge --allow-no-sex --filter-controls
#above command produces plink.frqx genotype count file for control samples

#output file from this code, 'assoc_1kgen_controls' can be used as an input instead of 'plink.assoc.adjusted' in other codes

#returns variants with p-value != NA into a list
def read_assoc(assoc_file):
    assoc_var = {}
    with open(assoc_file) as f:
        file=f.readlines()
    file = file[1:]
    for i in file:
        data = i.split()
        if data[8] != 'NA':
            assoc_var[data[1]] = data[4]
        else:
            pass
    return assoc_var


#adds control genotype counts to variable
def read_frqx(frqx_file, assoc_var):
    geno_count = {}
    with open(frqx_file) as f:
        file=f.readlines()
    file = file[1:]
    for i in file:
        data = i.split()
        if data[1] in assoc_var:
            case_min = assoc_var[data[1]]
            geno_count[data[1]] = data[:7] + [int(case_min)]
        else:
            pass
    return geno_count


def rest_api(rsid, pop_codes):
    pop_codes = pop_codes.upper().split(",")
    pop_codes = ["1000GENOMES:phase_3:"+i for i in pop_codes]
    count_dict = {}
    chrom = None
    alleles = None
    server = "https://rest.ensembl.org/variation/human/"
    ext = "?population_genotypes=1"
    url = server + rsid + ext
    headers={"User-Agent" : "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36", "Content-Type" : "application/json"}
    r = requests.get(url, headers=headers, verify=False)
    tries = 0
    while r.status_code != 200:
        tries += 1
        r = requests.get(url, headers=headers, verify=False)
        if tries > 100:
            break
    output = r.json()
    if ('population_genotypes' not in output or output['population_genotypes'] == []):
        pass
    else:
        alleles = output['mappings'][0]['allele_string'].split("/")
        chrom = output['mappings'][0]['seq_region_name']
        for val in output['population_genotypes']:
            if (val['population'] in pop_codes):
                geno = val['genotype']
                if geno not in count_dict:
                    count_dict[geno] = 0
                count_dict[geno] += val['count']
            else:
                pass
    return count_dict, chrom, alleles


def var_type(alleles):
    for ele in alleles:
        if (ele == "_" or len(ele)>1):
            var="indel"
            break
        else:
            var = "snp"
    return var


'''
Ensembl Rest_API and Plink allele naming can be different especially for indels.
Hence before adding the genotype counts they are compared and verified.
Multi-allelic sites are filtered.
For SNPs, the major alleles of both are compared, verified and added.
For indels with two alleles in the population, the long and short allele names are compared and verified based on the frequency.
For indels with single alleles, the major genotype count from rest-api is added to the plink major genotype.
Always combine genotype counts from same ethinicities and recheck the genotype counts of associated variants manually.
'''


def check_alleles(count_dict, value, chrom, alleles):
    var = var_type(alleles)
    plink_minor = value[2]
    plink_major = value[3]
    if len(plink_major) > len(plink_minor):
        input_major = "long_allele"
    else:
        input_major = "short_allele"
    zyg_rest = list(count_dict.keys())
    zyg_rest = [y for x in zyg_rest for y in x.split("|")]
    zyg_rest = list(set(zyg_rest))
    if len(zyg_rest) > 2:
        #variant is not bi-allelic
        res = "diff"
        rest_majgc, rest_mingc, het_count = 0, 0, 0
    else:
        if len(zyg_rest) == 1:
            #only single allele is present in the population
            if (zyg_rest[0] + "|" + zyg_rest[0]) not in count_dict:
                count_dict[zyg_rest[0] + "|" + zyg_rest[0]] = 0
            if zyg_rest[0] not in count_dict:
                count_dict[zyg_rest[0]] = 0
            #genotype count of single alleles are present in sex chromosomal variants extracted from ensembl
            #so half of single allele count is added to the genotype counts
            rest_majgc = count_dict[zyg_rest[0] + "|" + zyg_rest[0]] + ( 0.5 * count_dict[zyg_rest[0]] )
            rest_mingc = 0
            het_count = 0
            if var == "snp":
                rest_major = zyg_rest[0]
                if rest_major == plink_major:
                    res = "match"
                else:
                    res = "diff"  
            else:
                #variant is an indel
                #allele length of indels cannot be compared as it is a single allele
                #Hence the major genotype counts are added to the plink output
                res = "match"
        else:
            #both alleles are present in the populations
            for ele in zyg_rest:
                if (ele + "|" + ele) not in count_dict:
                    count_dict[ele + "|" + ele] = 0
                if ele not in count_dict:
                    count_dict[ele] = 0
            if (zyg_rest[0] + "|" + zyg_rest[1]) not in count_dict:
                count_dict[zyg_rest[0] + "|" + zyg_rest[1]] = 0
            if (zyg_rest[1] + "|" + zyg_rest[0]) not in count_dict:
                count_dict[zyg_rest[1] + "|" + zyg_rest[0]] = 0
            allele1 = zyg_rest[0]
            allele1_hcount = count_dict[zyg_rest[0] + "|" + zyg_rest[0]] + (0.5 * count_dict[zyg_rest[0]] )
            het_count = count_dict[zyg_rest[0] + "|" + zyg_rest[1]] +  count_dict[zyg_rest[1] + "|" + zyg_rest[0]]
            allele2 = zyg_rest[1]
            allele2_hcount = count_dict[zyg_rest[1] + "|" + zyg_rest[1]] + (0.5 * count_dict[zyg_rest[1]] )
            if allele1_hcount > allele2_hcount:
                rest_major = allele1
                rest_majgc = allele1_hcount
                rest_minor = allele2
                rest_mingc = allele2_hcount
            else:
                rest_major = allele2
                rest_majgc = allele2_hcount
                rest_minor = allele1
                rest_mingc = allele1_hcount
            if var == "snp":
                if rest_major == plink_major:
                    res = "match"
                else:
                    res = "diff"
            else:
                #variant is an indel
                if rest_major == "-":
                    rest_major = ""
                if rest_minor == "-":
                    rest_minor = ""
                if len(rest_major) > len(rest_minor):
                    rest_maj = "long_allele"
                else:
                    rest_maj = "short_allele" 
                if input_major == rest_maj:
                    res = "match"
                else:
                    res = "diff"
    return res, rest_majgc, rest_mingc, het_count


def read_file(var_file):
    var_data = []
    with open(var_file, 'r') as f:
        file = f.readlines()
    for i in file:
        single_var = i.split()[1]
        var_data.append(single_var)
    return var_data


def add_data(geno_count, pop_codes, var_data, var_file):
    new_list = []
    var_num = 0
    final_key = list(geno_count.keys())[-1]
    for key,value in geno_count.items():
        if key in var_data:
            continue
        #saves every 100 variant data into file
        if var_num >= 100:
            data_df = pd.DataFrame(new_list)
            print("Writing data of 100 variants to file...")
            data_df.to_csv(var_file, mode='a', sep=" ", index=False, header=False)
            new_list = []
            var_num = 0
        var_num += 1
        var_res = re.search(r'\brs\w+', key)
        if var_res == None:
            new_list.append(value)
            continue
        rsid = var_res.group()
        print("Extracting genotype counts of rsID:", rsid)
        count_dict, chrom, alleles = rest_api(rsid, pop_codes)
        if count_dict == {}:
            new_list.append(value)
            continue
        res, rest_majgc, rest_mingc, het_count = check_alleles(count_dict, value, chrom, alleles)
        if res == "match":
            var_data = value[:4] + [int(value[4])+rest_mingc , int(value[5])+het_count, int(value[6])+rest_majgc, value[7]]
            new_list.append(var_data)
        else:
            print("Allele validation did not match for:", rsid)
            new_list.append(value)
        if key == final_key:
            data_df = pd.DataFrame(new_list)
            print("Writing final variants to file...")
            data_df.to_csv(var_file, mode='a', sep=" ", index=False, header=False)
            #end of file


def hwe_test(hom_minor, het, hom_major):
    total = hom_minor + het + hom_major
    af_minor = (hom_minor + (0.5 * het)) / total
    af_major = (hom_major + (0.5 * het)) / total
    exp_gmin = (af_minor ** 2) * total
    exp_ghet = (2 * af_minor * af_major) * total
    exp_gmaj = (af_major ** 2) * total
    obs_list = [hom_minor, het, hom_major]
    exp_list = [exp_gmin, exp_ghet, exp_gmaj]
    #outputs chi-square test statistic and p-value
    cst, pval = chisquare(f_obs=obs_list,f_exp=exp_list)
    return cst, pval


def case_ctrl(case_min, case_max, ctrl_min, ctrl_max):
    row1_sum = case_min + ctrl_min
    row2_sum = case_max + ctrl_max
    col1_sum = case_min + case_max
    col2_sum = ctrl_min + ctrl_max
    total = row1_sum + row2_sum
    exp_camin = (row1_sum * col1_sum)/total
    exp_camax = (row2_sum * col1_sum)/total
    exp_ctmin = (row1_sum * col2_sum)/total
    exp_ctmax = (row2_sum * col2_sum)/total
    obs = [case_min, case_max, ctrl_min, ctrl_max]
    exp = [exp_camin, exp_camax, exp_ctmin, exp_ctmax]
    cst_assoc, pval_assoc = chisquare(f_obs=obs,f_exp=exp,ddof=2)
    #Haldane correction to avoid division by zero error.
    numer = (case_min + 0.5) / (ctrl_min + 0.5)
    denom = (case_max + 0.5) / (ctrl_max + 0.5)
    odds_ratio = numer / denom
    return cst_assoc, pval_assoc, odds_ratio

 
def assoc_test(var_file, thresh, case_num, assoc_num):
    with open(var_file, 'r') as f:
        new_file = f.readlines()
    assoc_list = []
    headers = ["CHR", "SNP", "UNADJ", "ODDS_RATIO", "MIN_ALLLE(A1)", "MAJOR_ALLELE(A2)", "HOM_A1", "HET", "HOM_A2", "CASE_A1", "CASE_A2"]
    for ele in new_file:
        hom_minor, het, hom_major = ele.split()[4:7]
        hom_minor, het, hom_major = float(hom_minor), float(het), float(hom_major)
        cst, pval = hwe_test(hom_minor, het, hom_major)
        if pval < thresh:
            pass
        else:
            case_min = float(ele.split()[7])
            case_max = (2 * case_num) - case_min
            ctrl_min = (hom_minor * 2) + het
            ctrl_max = (hom_major * 2) + het
            cst_assoc, pval_assoc, odds_ratio = case_ctrl(case_min, case_max, ctrl_min, ctrl_max)
            data = ele.split()[0:2] + [pval_assoc] + [odds_ratio] + ele.split()[2:7] + [case_min, case_max]
            assoc_list.append(data)
    assoc_df = pd.DataFrame(assoc_list, columns=headers)
    assoc_df["ADJ_PVAL"] = assoc_df["UNADJ"]
    if assoc_num == 0:
        assoc_df["ADJ_PVAL"] *= len(assoc_df)
    else:
        assoc_df["ADJ_PVAL"] *= assoc_num
    assoc_df = assoc_df.sort_values("UNADJ")
    return assoc_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add genotype counts to control data followed by association test')
    parser.add_argument("--assoc_num", type=int, help="number of variants in the input file, give value as 0 if entire plink.assoc file is used as input")
    parser.add_argument('--assoc_file', default = 'plink.assoc', help='plink association output file')
    parser.add_argument('--frqx_file', default = 'plink.frqx', help='plink genotype count report file')
    parser.add_argument('--thresh', type=int, default=1e-5, help='p-value threshold of Hardy-Weinberg equilibrium test for filtering out variants')
    parser.add_argument('--case_num', type=int, default=500, help='number of case samples')
    parser.add_argument('--out_file', default = 'assoc_1kgen_controls', help='output file with hwe filtered and adjusted associations')
    parser.add_argument('--pop_codes', default = 'ITU,STU', help='codes of 1000 genome phase-3 populations from which control genotype counts will be added, comma-separated')
    parser.add_argument("--var_file", default= 'var_file.txt', help="text file having variant data in which genotype counts are already added")
    args, unknown = parser.parse_known_args()
        
    print("Extracting variants from PLINK assoc output ...")
    assoc_var = read_assoc(args.assoc_file)
    
    print("Loading genotype counts from the control samples...")
    geno_count = read_frqx(args.frqx_file, assoc_var)
    
    #every batch of 100 variant data is saved to a file.
    #If the program crashes due to connectivity issues, running the code again will resume the search of genotype data using Rest-API after filtering already added variants.
    if os.path.isfile(args.var_file):
        print("Genotype counts of variants will be appended to the existing file...")
        var_data = read_file(args.var_file)
    else:
        print("Genotype counts of variants will be added to a new file...")
        var_data = []
    print("Adding control data from 1000 genome populations...")
    add_data(geno_count, args.pop_codes, var_data, args.var_file)
    
    print("Calculating HWE and Case/Control association results...")
    assoc_df = assoc_test(args.var_file, args.thresh, args.case_num, args.assoc_num)
    
    print("Writing updated association results to file...")
    assoc_df.to_csv(args.out_file, sep=" ", index=False)
    
