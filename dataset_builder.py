import os
import numpy as np
import pandas
import random

def enumerate_feature(feature_list):
    feature_options = set(feature_list)
    feature_option_num = {}
    counter = 0
    for c in feature_options:
        if c in feature_option_num:
            continue
        else:
            feature_option_num[c] = counter
            counter = counter + 1
    return feature_option_num

def create_dataset(disease, train_size=0, test_size=0, val_size=0):

    # ------- build positive dataset --------

    #location = input('Enter full path of datafile with associated SNPs: ')
    loc = "/Users/kavya/JHU/comp_bio/project/gwas-association-downloaded_2019-11-05-EFO_0001360-withChildTraits.tsv"
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

    context_option_num = enumerate_feature(contexts)

    positive_genes = set(file.MAPPED_GENE[indices].values)
    output_file_name = os.path.join(os.path.split(loc)[0], disease + '.tsv')
    with open(output_file_name, "w+") as out_file:
        header = "SNP\tCONTEXT\tINTERGENIC\tIS_RISK_FACTOR\n"
        out_file.write(header)
        for i in range(pos_train_size):
            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" + "1\n")

    # ------- build negative dataset -------

    loc = "/Users/kavya/JHU/comp_bio/project/gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv"
    file = pandas.read_csv(loc, sep='\t', lineterminator='\r', low_memory=False)
    total_snp_associations = len(file)
    
    indices = random.sample(range(total_snp_associations), neg_train_size * 2)
    snps = list(file.SNPS[indices].values)
    contexts = list(file.CONTEXT[indices].values)
    intergenic_info = list(file.INTERGENIC[indices].values)
    context_option_num = enumerate_feature(contexts)
    mapped_genes = list(file.MAPPED_GENE[indices].values)

    with open(output_file_name, "a") as out_file:
        counter = 0
        for i in range(neg_train_size * 2):
            #print(snps[i]) # possible formats for our snps is 's', 's; s; s;', 's x s'
            if (counter >= neg_train_size):
                break
            if (mapped_genes[i] in positive_genes):
                continue
            out_file.write(snps[i] + "\t" + str(context_option_num[contexts[i]]) + "\t" + str(intergenic_info[i]) + "\t" + "0\n")
            counter = counter + 1


def main():
    disease = "type2_diabetes"
    create_dataset(disease)
if __name__ == "__main__": main()