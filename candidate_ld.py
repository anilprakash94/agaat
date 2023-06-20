import pandas as pd
import re
import argparse

from candidate_gene import *
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
    header = ["CHR", "SNP", "BP", "A1", "C_A", "C_U", "A2", "CHISQ", "P", "OR", "GENE", "BONF"]
    for i in cls_var:
        line = i
        bonf = float(line[8]) * hap_blocks
        line += [bonf]
        adj_list.append(line)
    df = pd.DataFrame(adj_list, columns=header)
    sorted_df = df.sort_values("P")
    #No. of independent variants used for multiple-testing correction will be written at the end of the output file
    sorted_df.loc[len(sorted_df.index)] = ["Number of independent variants = " + str(hap_blocks)] + [""] * (len(header)-1)
    sorted_df.to_csv(out_file, index=False)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Candidate gene association analysis and multiple testing correcting using haplotype blocks')
    parser.add_argument('--block_file', default = 'plink.blocks', help='plink --blocks output file with haplotype blocks and variant IDs')
    parser.add_argument('--gene_list', default = 'gene_list.txt', help='input file having list of genes')
    parser.add_argument('--gene_bed', default = 'UCSC_canonical.bed', help='bed file with genomic coordinates of genes')
    parser.add_argument('--plink_assoc', default = 'plink.assoc', help='plink association output file')
    parser.add_argument('--out_file', default = 'plink.gene.ld.assoc.csv', help='output file with gene-list based associations')
    args, unknown = parser.parse_known_args()
    
    print("Extracting genes from text file...")
    cls_genes = get_genes(args.gene_list)
    
    print("Reading .bed file with gene coordinates")
    loc_dict = gene_coord(args.gene_bed, cls_genes)
    
    print("Retrieving plink variants that belong to the input list of genes...")
    cls_var = varto_gene(args.plink_assoc, loc_dict)
    
    print("Extracting variants and their corresponding ld blocks...")
    var_dict = block_var(args.block_file)
    
    print("Estimating number of independent markers using haplotype blocks...")
    hap_blocks = subset_var(cls_var, var_dict)
    
    print("Creating output csv file with candidate subset variants and their ld-based corrected p-value...")
    write_out(hap_blocks, cls_var, args.out_file)
    
