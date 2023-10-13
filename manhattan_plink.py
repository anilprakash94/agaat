import argparse, re
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text



def read_assoc(plink_assoc, var_num):
    with open(plink_assoc) as f:
        raw_p = f.readlines()
    var_list = []
    for i in range(1, len(raw_p)):
        var_info = raw_p[i].split()
        chrom, snp, pos = var_info[0:3]
        match_res = re.search(r'\brs\w+', snp)
        if match_res != None:
            snp = match_res.group()
        if chrom in ['23', '24', '25', '26']:
            continue
        pval = var_info[-2]
        if pval != 'NA':
            var_list.append([int(chrom), int(pos), snp, float(pval)])
    header = ["CHR", "POS", "SNP", "P-VALUE"]
    var_df = pd.DataFrame(var_list, columns=header)
    var_df.sort_values(by=['P-VALUE'],inplace=True)
    total_var = len(var_df)
    var_df = var_df.iloc[:var_num]
    var_df['-log10p'] = -np.log10(var_df['P-VALUE'])
    var_df.sort_values(by=['CHR','POS'],inplace=True)
    var_df['ind'] = range(len(var_df))
    return var_df, total_var


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plotting manhattan plot from PLINK output')
    parser.add_argument('--plink_assoc', default = 'plink.assoc', help='plink association file with raw p-values')
    parser.add_argument('--var_num', default = '2000', type=int, help='number of top variants to plot')  
    parser.add_argument('--plot_name', default = 'manhattan_plot.png', help='output plot name')
    args, unknown = parser.parse_known_args()
        
    print("Extracting PLINK output into dataframe...")
    var_df, total_var = read_assoc(args.plink_assoc, args.var_num)
    
    print("Plotting top significant variants...")
    
    plot = sns.relplot(data=var_df, x='ind', y='-log10p', aspect=4, hue='CHR', palette = 'bright', legend=False, height=6)
    bonf_pval = 0.05 / total_var
    threshold = -np.log10(bonf_pval)
    ax1 = plot.axes.flatten()[0]
    ax1.axhline(threshold, ls='--')
    chrom_df=var_df.groupby('CHR')['ind'].median()
    plot.ax.set_xlabel('Chromosome')
    plot.ax.set_xticks(chrom_df);
    plot.ax.set_xticklabels(chrom_df.index)
    bottom = ax1.get_ylim()[0]
    plot.set(ylim=(bottom,9))
    plot.fig.suptitle('Manhattan plot')
    #label significant SNPs
    sig_df = var_df[var_df['-log10p'] >= threshold]
    texts = []
    for var in range(len(sig_df)):
        x, y, label_txt = sig_df.iloc[var]['ind'], sig_df.iloc[var]['-log10p'], sig_df.iloc[var]['SNP']
        texts.append(ax1.annotate(label_txt,(x,y)))
    if texts != []:
        adjust_text(texts)
    plot.savefig(args.plot_name, dpi=300)
 


