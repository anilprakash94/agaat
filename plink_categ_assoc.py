import pandas as pd
import re
import argparse


def get_genes(gene_list):
    with open(gene_list) as f:
        gene_file=f.readlines()
    cls_genes = [i.strip() for i in gene_file]
    return cls_genes


def get_rsdata(dbsnp_common):
    with open(dbsnp_common) as vcf_file:
        vcf_data=vcf_file.readlines()
    #variant rows start from line 58
    vcf_data = vcf_data[57:]
    loc_data = {}
    for i in vcf_data:
        rsid = i.split()[2]
        res = re.search('GENEINFO\W(.+?)\:',i)
        if res:
            gene = res.group(1)
            loc_data[rsid] = gene
        else:
            pass
    return loc_data


def rsidto_gene(plink_adj, loc_data, cls_genes):
    with open(plink_adj) as f:
        file=f.readlines()
    file = file[1:]
    var_num = 0
    cls_var = []
    for i in file:
        var_num += 1
        print("Plink output variant no.:", var_num)
        mut = i.split()
        var = mut[1]
        unadj = mut[2]
        res = re.search(r'\brs\w+', var)
        if res:
            rsid = res.group()
            if rsid in loc_data:
                gene = loc_data[rsid]
                for name in gene.split("-"):
                    if name in cls_genes:
                        cls_var.append(mut[:3]+[gene])
                        break
            else:
                pass     
        else:
            pass
    return cls_var
 

def create_csv(cls_var, out_file):
    adj_list = []
    header = ["chr", "rsID", "unadj", "gene", "bonf"]
    adj_list.append(header)
    for ele in cls_var:
        bonf = float(ele[2]) * len(cls_var)
        line = ele + [bonf]
        adj_list.append(line)
    df = pd.DataFrame(adj_list)
    adj_df = df.rename(columns=df.iloc[0]).loc[1:]
    adj_df.to_csv(out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Multiple testing correction using genetic variants belonging to a category of genes')
    parser.add_argument('--gene_list', default = 'gene_list.txt', help='input file having list of genes')
    parser.add_argument('--dbsnp_common', default = 'common_all_20180418.vcf', help='dbsnp vcf file with common variants')
    parser.add_argument('--plink_adj', default = 'plink.assoc.adjusted', help='plink file with adjusted associations')
    parser.add_argument('--out_file', default = 'plink.category.adjusted.csv', help='output file with gene-list based adjusted associations')
    args, unknown = parser.parse_known_args()

    print("Extracting gene list from text file...")
    cls_genes = get_genes(args.gene_list)
    
    print("Extracting rsIDs of common variants and their corresponding genes from dbSNP common_all_20180418.vcf dataset...")
    loc_data = get_rsdata(args.dbsnp_common)
    
    print("Retrieving plink variants that belong to the input list of genes...")
    cls_var = rsidto_gene(args.plink_adj, loc_data, cls_genes)
    
    print("Creating output csv file with plink variants, genes and their corrected p-value...")
    create_csv(cls_var, args.out_file)

