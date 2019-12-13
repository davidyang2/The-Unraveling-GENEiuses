import math
import matplotlib.pyplot as plt
import numpy as np


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


def get_feature_counts_for_plotting(feature_list, feature_option_num={}):
    counts = {}
    for f in feature_list:
        if feature_option_num == {}:
            f_num = f
        else:
            f_num = feature_option_num[f]
        if f_num in counts:
            counts[f_num] = counts[f_num] + 1
        else:
            counts[f_num] = 0
    return counts


def plot_feature(positive_counts, negative_counts, feature_type):
    width = 0.35
    fig, ax = plt.subplots()
    lower_range = min(min(positive_counts), min(negative_counts))
    upper_range = max(max(positive_counts), max(negative_counts))
    ind = max(len(positive_counts), len(negative_counts))
    #print(ind)

    if feature_type == "GO FUNCTION":
        ax.set_xlim(-1, ind+1)
    else:
        ax.set_xlim(lower_range - 1, upper_range + 1)

    # pp = PdfPages(feature_type + "_histograms.pdf")
    # ax = plt.gca()
    # ax.set_xlim(0, 18)
    # pos_context_fig = plt.figure(1)
    ind = np.fromiter(positive_counts.keys(), dtype=float)
    pos_context_fig = ax.bar(ind, positive_counts.values(), width, color='g')
    # pp.savefig(pos_context_fig, dpi = 300, transparent=True)

    # neg_context_fig = plt.figure(2)
    ind = np.fromiter(negative_counts.keys(), dtype=float)
    neg_context_fig = ax.bar(ind + width, negative_counts.values(), width, color='r')
    # pp.savefig(neg_context_fig, dpi=300, transparent=True)
    # pp.close()

    ax.set_title(feature_type)
    ax.set_ylabel('Count')
    ax.set_xlabel('Type')

    fig.savefig(feature_type + "_histograms.pdf")


def error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance,
                                               downstream_gene_distance, reported_genes, mapped_genes):
    if math.isnan(intergenic_info[i]):
        return True
    try:
        chromosome_ids[i] = float(chromosome_ids[i])
        chr_pos[i] = float(chr_pos[i])
        upstream_gene_distance[i] = float(upstream_gene_distance[i])
        downstream_gene_distance[i] = float(downstream_gene_distance[i])
    except ValueError:
        return True

    if (chromosome_ids[i] == "") or (chr_pos[i] == "") or (upstream_gene_distance[i] == ""):
        return True

    if math.isnan(chromosome_ids[i]) or math.isnan(chr_pos[i]):
        return True

    if reported_genes[i] == "NR" or reported_genes[i] == "":
        return True

    if mapped_genes[i] == "":
        return True

    return False


def enumerate_gwas_features(contexts, intergenic_info, context_option_num={}, intergenic_option_num={}):
    context_option_num = enumerate_feature(contexts, context_option_num)
    intergenic_option_num = enumerate_feature(intergenic_info, intergenic_option_num)
    return context_option_num, intergenic_option_num

def enumerate_a_feature(feature, feature_option_num={}):
    feature_option_num = enumerate_feature(contexts, context_option_num)
    return feature_option_num


def get_gwas_features(file, indices):
    contexts = list(file.CONTEXT[indices].values)
    intergenic_info = list(file.INTERGENIC[indices].values)

    # new features
    chromosome_ids = list(file.CHR_ID[indices].values)
    chr_pos = list(file.CHR_POS[indices].values)
    upstream_gene_distance = list(file.UPSTREAM_GENE_DISTANCE[indices].values)
    downstream_gene_distance = list(file.DOWNSTREAM_GENE_DISTANCE[indices].values)

    # For SNPs that are not intergenic, GWAS lists their upstream and downstream gene distance as -1
    for i in range(len(upstream_gene_distance)):
        if math.isnan(upstream_gene_distance[i]):
            upstream_gene_distance[i] = -1

    for j in range(len(downstream_gene_distance)):
        if math.isnan(downstream_gene_distance[j]):
            downstream_gene_distance[j] = -1

    return contexts, intergenic_info, chromosome_ids, chr_pos, upstream_gene_distance, downstream_gene_distance
