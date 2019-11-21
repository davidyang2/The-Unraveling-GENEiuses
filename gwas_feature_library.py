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

def error_raised_in_gwas_features_for_curr_snp(i, intergenic_info, chromosome_ids, chr_pos, \
    upstream_gene_distance, downstream_gene_distance):

    if (math.isnan(intergenic_info[i])):
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

    if (math.isnan(chromosome_ids[i]) or math.isnan(chr_pos[i]) or math.isnan(upstream_gene_distance[i]) \
        or math.isnan(downstream_gene_distance[i])):
        return True

    return False

def enumerate_gwas_features(contexts, intergenic_info, \
    context_option_num={}, intergenic_option_num={}):
    context_option_num = enumerate_feature(contexts, context_option_num)
    intergenic_option_num = enumerate_feature(intergenic_info, intergenic_option_num)
    return context_option_num, intergenic_option_num

def get_gwas_features(file, indices):
    contexts = list(file.CONTEXT[indices].values)
    intergenic_info = list(file.INTERGENIC[indices].values)

    # new features
    chromosome_ids = list(file.CHR_ID[indices].values)
    chr_pos = list(file.CHR_POS[indices].values)
    upstream_gene_distance = list(file.UPSTREAM_GENE_DISTANCE[indices].values)
    downstream_gene_distance = list(file.DOWNSTREAM_GENE_DISTANCE[indices].values)
    return contexts, intergenic_info, chromosome_ids, chr_pos, \
        upstream_gene_distance, downstream_gene_distance