import os
import numpy as np
import pandas
import random
import math
import sys

import matplotlib.pyplot as plt

import gwas_feature_library

def build_positive_dataset(loc, output_file_name, train_size=0):
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
    positive_genes = set(file.MAPPED_GENE[indices].values)

    contexts, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, \
        downstream_gene_distance = gwas_feature_library.get_gwas_features(file, indices)

    context_option_num, intergenic_option_num = gwas_feature_library.enumerate_gwas_features(contexts, \
        intergenic_info)

    positive_context_counts = gwas_feature_library.get_feature_counts_for_plotting(contexts, \
        context_option_num)
    positive_intergenic_counts = gwas_feature_library.get_feature_counts_for_plotting(intergenic_info, \
        intergenic_option_num)

    with open(output_file_name, "w+") as out_file:
        header = "SNP\tCONTEXT\tINTERGENIC\tCHR_ID\tCHR_POS\tUP_DIST\tDOWN_DIST\t"

        # Example code:
        # header += kavya_features.add_to_header()
        # header += david_features.add_to_header()
        # header += richard_features.add_to_header()

        header += "IS_RISK_FACTOR\n"
        
        out_file.write(header)
        for i in range(pos_train_size):

            if (gwas_feature_library.error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, chromosome_ids, \
                chr_pos, upstream_gene_distance, downstream_gene_distance)):
                continue

            try:
                # Example code
                features_to_add = ""
                # features_to_add += kavya_features.add_feature_for_curr_snp(snps[i])
                # features_to_add += david_features.add_feature_for_curr_snp(snps[i])
                # features_to_add += richard_features.add_feature_for_curr_snp(snps[i])
            except ValueError:
                continue

            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" \
                + str(chromosome_ids[i]) + "\t"  + str(chr_pos[i]) + "\t" \
                    + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t")

            # Example code:
            # out_file.write(features_to_add)

            out_file.write("1\n")

    return pos_train_size, positive_genes, context_option_num, intergenic_option_num, \
        positive_context_counts, positive_intergenic_counts


def build_negative_dataset(loc, output_file_name, neg_train_size, positive_genes, \
    context_option_num, intergenic_option_num, positive_context_counts, positive_intergenic_counts):

    file = pandas.read_csv(loc, sep='\t', lineterminator='\r', low_memory=False)
    total_snp_associations = len(file)
    
    indices = random.sample(range(total_snp_associations), neg_train_size)
    snps = list(file.SNPS[indices].values)
    mapped_genes = list(file.MAPPED_GENE[indices].values)

    contexts, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, \
        downstream_gene_distance = gwas_feature_library.get_gwas_features(file, indices)

    context_option_num, intergenic_option_num = gwas_feature_library.enumerate_gwas_features(contexts, \
        intergenic_info, context_option_num, intergenic_option_num)

    negative_intergenic_counts = gwas_feature_library.get_feature_counts_for_plotting(intergenic_info, \
        intergenic_option_num)
    negative_context_counts = gwas_feature_library.get_feature_counts_for_plotting(contexts, \
        context_option_num)

    gwas_feature_library.plot_feature(positive_context_counts, negative_context_counts, "Context")
    gwas_feature_library.plot_feature(positive_intergenic_counts, negative_intergenic_counts, "Intergenic")

    with open(output_file_name, "a") as out_file:
        counter = 0
        for i in range(neg_train_size):
            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            if (counter >= neg_train_size):
                break
            if (mapped_genes[i] in positive_genes):
                continue

            if (gwas_feature_library.error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, \
                chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance)):
                continue

            try:
                # Example code
                features_to_add = ""
                # features_to_add += kavya_features.add_feature_for_curr_snp(snps[i])
                # features_to_add += david_features.add_feature_for_curr_snp(snps[i])
                # features_to_add += richard_features.add_feature_for_curr_snp(snps[i])
            except ValueError:
                continue

            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" \
                + str(chromosome_ids[i]) + "\t"  + str(chr_pos[i]) + "\t" \
                    + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t")

            # Example code:
            # out_file.write(features_to_add)

            out_file.write("0\n")

            counter = counter + 1

def create_dataset(disease, train_size=0, test_size=0, val_size=0):
    loc = "/Users/kavya/JHU/comp_bio/project/gwas-diabetes2-positive.tsv"
    out_file_name = os.path.join(os.path.split(loc)[0], disease + '.tsv')

    pos_train_size, pos_genes, context_option_num, intergenic_option_num, positive_context_counts, \
        positive_intergenic_counts = build_positive_dataset(loc, out_file_name)

    neg_train_size = pos_train_size
    loc = "/Users/kavya/JHU/comp_bio/project/gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv"
    build_negative_dataset(loc, out_file_name, neg_train_size, pos_genes, context_option_num, \
        intergenic_option_num, positive_context_counts, positive_intergenic_counts)


def main():
    disease = "diabetes"
    create_dataset(disease)
if __name__ == "__main__": main()