import argparse
import pandas as pd


def block_var(block_file):
    with open(block_file) as f:
        file = f.readlines()
    var_dict = {}
    for i in range(len(file)):
        var_list = file[i][1:].strip().replace(',',' ').split()
        for ele in var_list:
            var_dict[ele] = i
    return var_dict


def assoc_var(assoc_adj, var_dict):
    block_num = []
    ind_var = 0
    with open(assoc_adj) as f:
        assoc_file = f.readlines()
    assoc_file = assoc_file[1:]
    for i in assoc_file:
        varnt = i.split()[1]
        if varnt in var_dict:
            block_num.append(var_dict[varnt])
        else:
            ind_var += 1
    hap_blocks = len(set(block_num))
    hap_blocks += ind_var
    return hap_blocks, assoc_file


def write_out(hap_blocks, assoc_file, out_file):
    adj_list = []
    header = ["chr", "rsID", "unadj", "bonf"]
    for i in assoc_file:
        line = i.split()[:3]
        bonf = float(line[2]) * hap_blocks
        line += [bonf]
        adj_list.append(line)
    df = pd.DataFrame(adj_list, columns=header)
    df.to_csv(out_file, index=False)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Multiple testing correction using haplotype blocks')
    parser.add_argument('--block_file', default = 'plink.blocks', help='plink --blocks output file with haplotype blocks and variant IDs')
    parser.add_argument('--assoc_adj', default = 'plink.assoc.adjusted', help='plink file with adjusted associations')
    parser.add_argument('--out_file', default = 'plink.ld.adjusted.csv', help='output file with ld-block based adjusted associations')
    args, unknown = parser.parse_known_args()
        
    print("Extracting variants and their corresponding ld blocks...")
    var_dict = block_var(args.block_file)
    
    print("Estimating number of independent markers using haplotype blocks...")
    hap_blocks, assoc_file = assoc_var(args.assoc_adj, var_dict)
    
    print("Creating output csv file with associated variants and their corrected p-value...")
    write_out(hap_blocks, assoc_file, args.out_file)

