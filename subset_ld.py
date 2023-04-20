import pandas as pd
import re
import argparse

from plink_categ_assoc import *
from ld_blocks import *


def subset_var(cls_var, var_dict):
    block_num = []
    ind_var = 0
    for i in cls_var:
        varnt = i[1]
        if varnt in var_dict:
            block_num.append(var_dict[varnt])
        else:
            ind_var += 1
    hap_blocks = len(set(block_num))
    hap_blocks += ind_var
    return hap_blocks


def write_out(hap_blocks, cls_var, out_file):
    adj_list = []
    header = ["chr", "rsID", "unadj", "gene", "bonf"]
    for i in cls_var:
        line = i
        bonf = float(line[2]) * hap_blocks
        line += [bonf]
        adj_list.append(line)    
    #No. of independent variants used for multiple-testing correction will be written at the end of the output file
    adj_list.append(["Number of independent variants = " + str(hap_blocks), "", "", "", ""])
    df = pd.DataFrame(adj_list, columns=header)
    df.to_csv(out_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Candidate gene association analysis and multiple testing correcting using haplotype blocks')
    parser.add_argument('--block_file', default = 'plink.blocks', help='plink --blocks output file with haplotype blocks and variant IDs')
    parser.add_argument('--gene_list', default = 'gene_list.txt', help='input file having list of genes')
    parser.add_argument('--dbsnp_common', default = 'common_all_20180418.vcf', help='dbsnp vcf file with common variants')
    parser.add_argument('--plink_adj', default = 'plink.assoc.adjusted', help='plink file with adjusted associations')
    parser.add_argument('--out_file', default = 'plink.subset.ld.csv', help='output file with ld-block based adjusted associations')
    args, unknown = parser.parse_known_args()

    print("Extracting gene list from text file...")
    cls_genes = get_genes(args.gene_list)
    
    print("Extracting rsIDs of common variants and their corresponding genes from dbSNP common_all_20180418.vcf dataset...")
    loc_data = get_rsdata(args.dbsnp_common)
    
    print("Retrieving output variants from plink that belong to the input list of genes...")
    cls_var = rsidto_gene(args.plink_adj, loc_data, cls_genes)
    
    print("Extracting variants and their corresponding ld blocks...")
    var_dict = block_var(args.block_file)
    
    print("Estimating number of independent markers using haplotype blocks...")
    hap_blocks = subset_var(cls_var, var_dict)
    
    print("Creating output csv file with candidate subset variants and their ld-based corrected p-value...")
    write_out(hap_blocks, cls_var, args.out_file)
    
