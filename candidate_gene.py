import pandas as pd
import argparse


def get_genes(gene_list):
    with open(gene_list) as f:
        gene_file=f.readlines()
    cls_genes = {i.strip().upper() for i in gene_file}
    return cls_genes


#reads .bed file with gene coordinates

def gene_coord(gene_bed, cls_genes):
    with open(gene_bed) as f:
        file = f.readlines()
    file = file[1:]
    loc_dict = {}
    for i in file:
        ele = i.split()
        gene_name = ele[4].upper()
        if gene_name not in cls_genes:
            continue
        if ele[0] not in loc_dict:
            loc_dict[ele[0]] = []
        loc_dict[ele[0]].append(ele[1]+":"+ele[2]+":"+gene_name)
    return loc_dict


def varto_gene(plink_assoc, loc_dict):
    with open(plink_assoc) as f:
        file=f.readlines()
    file = file[1:]
    plink_annot = {"23" : "X", "24" : "Y", "25" : "X", "26" : "MT"}
    var_num = 0
    cls_var = []
    for i in file:
        var_num += 1
        print("Plink output variant no.:", var_num)
        mut = i.split()
        if mut[8] == "NA":
            continue
        chr_num = mut[0]
        if chr_num in plink_annot:
            chr_num = plink_annot[chr_num]
        chrom = "chr" + chr_num
        if chrom not in loc_dict:
            continue
        coord = float(mut[2])
        for ele in loc_dict[chrom]:
            gene_coord = ele.split(":")
            if coord >= float(gene_coord[0]) and coord <= float(gene_coord[1]):
                gene_name = gene_coord[2]
                cls_var.append(mut + [gene_name])
                break
            else:
                pass
    return cls_var
 

def create_csv(cls_var, out_file):
    adj_list = []
    header = ["CHR", "SNP", "BP", "A1", "C_A", "C_U", "A2", "CHISQ", "P", "OR", "GENE", "BONF"]
    for ele in cls_var:
        bonf = float(ele[8]) * len(cls_var)
        line = ele + [bonf]
        adj_list.append(line)
    df = pd.DataFrame(adj_list, columns=header)
    sorted_df = df.sort_values("P")
    sorted_df.to_csv(out_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Multiple testing correction using genetic variants belonging to a category of genes')
    parser.add_argument('--gene_list', default = 'gene_list.txt', help='input file having list of genes')
    parser.add_argument('--gene_bed', default = 'UCSC_canonical.bed', help='bed file with genomic coordinates of genes')
    parser.add_argument('--plink_assoc', default = 'plink.assoc', help='plink association output file')
    parser.add_argument('--out_file', default = 'plink.gene.assoc.csv', help='output file with gene-list based associations')
    args, unknown = parser.parse_known_args()

    print("Extracting genes from text file into a set...")
    cls_genes = get_genes(args.gene_list)
    
    print("Reading .bed file with gene coordinates")
    loc_dict = gene_coord(args.gene_bed, cls_genes)
    
    print("Retrieving plink variants that belong to the input list of genes...")
    cls_var = varto_gene(args.plink_assoc, loc_dict)
    
    print("Creating output csv file with plink variants, genes and their corrected p-value...")
    create_csv(cls_var, args.out_file)

