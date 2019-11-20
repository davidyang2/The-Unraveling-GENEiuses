import os
import numpy as np
import pandas
import random
import math
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def enumerate_feature(feature_list, feature_option_num):
    feature_options = set(feature_list)
    counter = 0
    for c in feature_options:
        if c in feature_option_num:
            continue
        else:
            feature_option_num[c] = counter
            counter = counter + 1
    return feature_option_num

def get_feature_counts_for_plotting(feature_list, feature_option_num):
    counts = {}
    for f in feature_list:
        f_num = feature_option_num[f]
        if f_num in counts:
            counts[f_num] = counts[f_num] + 1
        else:
            counts[f_num] = 0
    return counts

def plot_feature(positive_counts, negative_counts, feature_type):
    width = 0.35
    fig, ax = plt.subplots()
    ind = max(len(positive_counts), len(negative_counts))
    ax.set_xlim(-1, 2)

    #pp = PdfPages(feature_type + "_histograms.pdf")
    #ax = plt.gca()
    #ax.set_xlim(0, 18)
    #pos_context_fig = plt.figure(1)
    ind = np.fromiter(positive_counts.keys(), dtype=float)
    pos_context_fig = ax.bar(ind, positive_counts.values(), width, color='g')
    #pp.savefig(pos_context_fig, dpi = 300, transparent=True)

    #neg_context_fig = plt.figure(2)
    ind = np.fromiter(negative_counts.keys(), dtype=float)
    neg_context_fig = ax.bar(ind + width, negative_counts.values(), width, color='r')
    #pp.savefig(neg_context_fig, dpi=300, transparent=True)
    #pp.close()

    ax.set_title(feature_type)
    ax.set_ylabel('Count')
    ax.set_xlabel('Type')

    fig.savefig(feature_type + "_histograms.pdf")


def create_dataset(disease, train_size=0, test_size=0, val_size=0):

    # ------- build positive dataset --------

    #location = input('Enter full path of datafile with associated SNPs: ')
    #loc = sys.argv[1]
    loc = "/Users/kavya/JHU/comp_bio/project/gwas-diabetes2-positive.tsv"
    file = pandas.read_csv(loc, sep='\t', lineterminator='\r')
    total_snp_associations = len(file)

    if train_size == 0:
        train_size = total_snp_associations
    pos_train_size = train_size // 2
    neg_train_size = train_size // 2

    # select snps to be used in training dataset
    random.seed(100)
    indices = random.sample(range(total_snp_associations), pos_train_size)
    snps = list(file.SNPS[indices].values)
    contexts = list(file.CONTEXT[indices].values)
    intergenic_info = list(file.INTERGENIC[indices].values)

    # new features
    chromosome_ids = list(file.CHR_ID[indices].values)
    chr_pos = list(file.CHR_POS[indices].values)
    upstream_gene_distance = list(file.UPSTREAM_GENE_DISTANCE[indices].values)
    downstream_gene_distance = list(file.DOWNSTREAM_GENE_DISTANCE[indices].values)

    context_option_num = enumerate_feature(contexts, {})
    positive_context_counts = get_feature_counts_for_plotting(contexts, context_option_num)

    intergenic_option_num = enumerate_feature(intergenic_info, {})
    positive_intergenic_counts = get_feature_counts_for_plotting(intergenic_info, intergenic_option_num)

    #chromosome_option_num = enumerate_feature(chromosome_ids, {})

    #print("POSITIVE")
    #for i in context_option_num:
        #print(str(i) + " " + str(context_option_num[i]))

    positive_genes = set(file.MAPPED_GENE[indices].values)
    output_file_name = os.path.join(os.path.split(loc)[0], disease + '.tsv')
    with open(output_file_name, "w+") as out_file:
        header = "SNP\tCONTEXT\tINTERGENIC\tCHR_ID\tCHR_POS\tUP_DIST\tDOWN_DIST\tIS_RISK_FACTOR\n"
        out_file.write(header)
        for i in range(pos_train_size):
            if (math.isnan(intergenic_info[i])):
                continue
            try:
                chromosome_ids[i] = float(chromosome_ids[i])
                chr_pos[i] = float(chr_pos[i])
                upstream_gene_distance[i] = float(upstream_gene_distance[i])
                downstream_gene_distance[i] = float(downstream_gene_distance[i])
            except ValueError:
                continue

            if (chromosome_ids[i] == "") or (chr_pos[i] == "") or (upstream_gene_distance[i] == ""):
                continue

            if (math.isnan(chromosome_ids[i]) or math.isnan(chr_pos[i]) or math.isnan(upstream_gene_distance[i]) \
                or math.isnan(downstream_gene_distance[i])):
                continue

            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" \
                + str(chromosome_ids[i]) + "\t"  + str(chr_pos[i]) + "\t" \
                    + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t" + "1\n")

    # ------- build negative dataset -------

    #loc = sys.argv[2]
    loc = "/Users/kavya/JHU/comp_bio/project/gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv"
    file = pandas.read_csv(loc, sep='\t', lineterminator='\r', low_memory=False)
    total_snp_associations = len(file)
    
    indices = random.sample(range(total_snp_associations), neg_train_size)
    snps = list(file.SNPS[indices].values)
    contexts = list(file.CONTEXT[indices].values)
    intergenic_info = list(file.INTERGENIC[indices].values)
    mapped_genes = list(file.MAPPED_GENE[indices].values)
    chromosome_ids = list(file.CHR_ID[indices].values)

    # new features
    chromosome_ids = list(file.CHR_ID[indices].values)
    chr_pos = list(file.CHR_POS[indices].values)
    upstream_gene_distance = list(file.UPSTREAM_GENE_DISTANCE[indices].values)
    downstream_gene_distance = list(file.DOWNSTREAM_GENE_DISTANCE[indices].values)

    context_option_num = enumerate_feature(contexts, context_option_num)
    intergenic_option_num = enumerate_feature(intergenic_info, intergenic_option_num)
    #chromosome_option_num = enumerate_feature(chromosome_ids, chromosome_option_num)

    negative_intergenic_counts = get_feature_counts_for_plotting(intergenic_info, intergenic_option_num)
    negative_context_counts = get_feature_counts_for_plotting(contexts, context_option_num)

    plot_feature(positive_context_counts, negative_context_counts, "Context")
    plot_feature(positive_intergenic_counts, negative_intergenic_counts, "Intergenic")

    #print("\n\n\nNEGATIVE")
    #for i in context_option_num:
        #print(str(i) + " " + str(context_option_num[i]))

    with open(output_file_name, "a") as out_file:
        counter = 0
        for i in range(neg_train_size):
            if (math.isnan(intergenic_info[i])):
                continue
            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            if (counter >= neg_train_size):
                break
            if (mapped_genes[i] in positive_genes):
                continue

            try:
                chromosome_ids[i] = float(chromosome_ids[i])
                chr_pos[i] = float(chr_pos[i])
                upstream_gene_distance[i] = float(upstream_gene_distance[i])
                downstream_gene_distance[i] = float(downstream_gene_distance[i])
            except ValueError:
                continue

            if (math.isnan(chromosome_ids[i]) or math.isnan(chr_pos[i]) or math.isnan(upstream_gene_distance[i]) \
                or math.isnan(downstream_gene_distance[i])):
                continue

            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" \
                + str(chromosome_ids[i]) + "\t"  + str(chr_pos[i]) + "\t" \
                    + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t" + "0\n")

            counter = counter + 1


def main():
    disease = "diabetes"
    create_dataset(disease)
if __name__ == "__main__": main()