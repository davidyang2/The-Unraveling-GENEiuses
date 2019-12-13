import os
import numpy as np
import pandas
import random
import math
import sys

import matplotlib.pyplot as plt

import gene_tissue_expression
import gwas_feature_library
import gene_ontology_features


def select_indices(pos_loc, neg_loc):
    file = pandas.read_csv(pos_loc, sep='\t', lineterminator='\r')
    positive_snp_associations = len(file)

    train_size = positive_snp_associations
    pos_train_size = train_size  # Size of positive data set is the same size as the length of the positive dataset file
    neg_train_size = train_size  # Size of negative data set is the same size as the length of the positive dataset file

    # select snps to be used in training dataset
    random.seed(100)
    pos_indices = random.sample(range(positive_snp_associations), pos_train_size)  # Randomly selects rows in the positive dataset file

    file = pandas.read_csv(neg_loc, sep='\t', lineterminator='\r', low_memory=False)
    negative_snp_associations = len(file)

    neg_indices = random.sample(range(negative_snp_associations), neg_train_size)  # Randomly selects rows in the negative dataset file

    return pos_indices, neg_indices


def build_positive_dataset(disease, loc, indices, output_file_name, train_size=0):
    file = pandas.read_csv(loc, sep='\t', lineterminator='\r')
    pos_train_size = len(indices)

    snps = list(file.SNPS[indices].values)
    positive_genes = set(file.MAPPED_GENE[indices].values)

    contexts, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance = gwas_feature_library.get_gwas_features(file, indices)

    context_option_num, intergenic_option_num = gwas_feature_library.enumerate_gwas_features(contexts, intergenic_info)

    positive_context_counts = gwas_feature_library.get_feature_counts_for_plotting(contexts, context_option_num)
    positive_intergenic_counts = gwas_feature_library.get_feature_counts_for_plotting(intergenic_info, intergenic_option_num)

    df = pandas.DataFrame(file)
    reported_genes = df[df.columns[13]]
    reported_genes = reported_genes[indices].values
    mapped_genes = df["MAPPED_GENE"].values

    with open(output_file_name, "w+") as out_file:
        header = "SNP\tCONTEXT\tINTERGENIC\tCHR_ID\tCHR_POS\tUP_DIST\tDOWN_DIST\t"

        # Example code:
        header += gene_ontology_features.add_to_header()
        header += gene_tissue_expression.add_to_header(disease)
        # header += richard_features.add_to_header()

        header += "IS_RISK_FACTOR\n"

        out_file.write(header)
        for i in range(pos_train_size):

            if gwas_feature_library.error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance, reported_genes, mapped_genes):
                continue
            try:
                # Example code
                features_to_add = ""
                # features_to_add += gene_ontology_features.add_feature_for_curr_snp(disease, reported_genes[i])
                if disease == "cardiovascular_disease":
                    features_to_add += gene_tissue_expression.add_tissue_expression_for_curr_snp(disease, mapped_genes[i])
                else:
                    features_to_add += gene_tissue_expression.add_tissue_expression_for_curr_snp(disease, reported_genes[i])
                print(features_to_add)
            except ValueError:
                continue

            # print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" + str(chromosome_ids[i]) + "\t" + str(chr_pos[i]) + "\t" + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t")

            # Example code:
            print("posWrite")
            out_file.write(features_to_add)

            out_file.write("1\n")

    return positive_genes, context_option_num, intergenic_option_num, positive_context_counts, positive_intergenic_counts


def build_negative_dataset(disease, loc, indices, output_file_name, positive_genes, context_option_num,
                           intergenic_option_num, positive_context_counts, positive_intergenic_counts):
    file = pandas.read_csv(loc, sep='\t', lineterminator='\r', low_memory=False)
    neg_train_size = len(indices)

    snps = list(file.SNPS[indices].values)
    # mapped_genes = list(file.MAPPED_GENE[indices].values)

    contexts, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance = gwas_feature_library.get_gwas_features(file, indices)

    context_option_num, intergenic_option_num = gwas_feature_library.enumerate_gwas_features(contexts, intergenic_info, context_option_num, intergenic_option_num)

    negative_intergenic_counts = gwas_feature_library.get_feature_counts_for_plotting(intergenic_info, intergenic_option_num)
    negative_context_counts = gwas_feature_library.get_feature_counts_for_plotting(contexts, context_option_num)

    gwas_feature_library.plot_feature(positive_context_counts, negative_context_counts, "Context")
    gwas_feature_library.plot_feature(positive_intergenic_counts, negative_intergenic_counts, "Intergenic")

    df = pandas.DataFrame(file)
    reported_genes = df["REPORTED GENE(S)"] # df[df.columns[13]]
    reported_genes = reported_genes[indices].values
    mapped_genes = df["MAPPED_GENE"].values

    with open(output_file_name, "a") as out_file:
        counter = 0
        for i in range(neg_train_size):
            # print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            if counter >= neg_train_size:
                break
            if mapped_genes[i] in positive_genes:
                continue

            if gwas_feature_library.error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance, reported_genes, mapped_genes):
                continue
            try:
                # Example code
                features_to_add = ""
                features_to_add += gene_ontology_features.add_feature_for_curr_snp(disease, reported_genes[i])
                if disease == "cardiovascular_disease":
                    features_to_add += gene_tissue_expression.add_tissue_expression_for_curr_snp(disease, mapped_genes[i])
                else:
                    features_to_add += gene_tissue_expression.add_tissue_expression_for_curr_snp(disease, reported_genes[i])
            except ValueError:
                continue

            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t"
                           + str(intergenic_info[i]) + "\t" + str(chromosome_ids[i]) + "\t" + str(chr_pos[i])
                           + "\t" + str(upstream_gene_distance[i]) + "\t" + str(downstream_gene_distance[i]) + "\t")

            # Example code:
            print("negWrite")
            out_file.write(features_to_add)

            out_file.write("0\n")

            counter = counter + 1


def create_dataset(disease, train_size=0, test_size=0, val_size=0):
    pos_loc = sys.argv[1]  # C:\Users\david\Documents\Computational Biomedical Research\gwas-association-downloaded_2019-11-06-EFO_0000305-withChildTraits.tsv
    neg_loc = sys.argv[2]  # C:\Users\david\Documents\Computational Biomedical Research\gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv

    pos_indices, neg_indices = select_indices(pos_loc, neg_loc)  # Gets rows for genes of interest in positive and negative dataset

    gene_ontology_features.make_gene_list(disease, pos_loc, neg_loc, pos_indices, neg_indices)  # Gene list sent to a txt file for Gene Ontology

    out_file_name = os.path.join(os.path.split(pos_loc)[0], disease + '.tsv')

    pos_genes, context_option_num, intergenic_option_num, positive_context_counts, positive_intergenic_counts = build_positive_dataset(disease, pos_loc, pos_indices, out_file_name)

    build_negative_dataset(disease, neg_loc, neg_indices, out_file_name, pos_genes, context_option_num, intergenic_option_num, positive_context_counts, positive_intergenic_counts)


def main():
    disease = "cardiovascular_disease"  # Change this depending on what disease you are analyzing
    create_dataset(disease)


if __name__ == "__main__": main()
