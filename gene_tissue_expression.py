import pandas as pd
import os
import sys


def add_to_header():
    return "TISSUE_EXPRESSION LEVEL\t"

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

def add_tissue_expression_for_curr_snp(disease, curr_gene):
    expression_dict = {"Not detected": "0", "Low": "1", "Medium": "2", "High": "3"}

    curr_gene = cleanup_single_gene(str(curr_gene))
    expression_report_loc = "/Users/david/Documents/Computational Biomedical Research/normal_tissue.tsv/normal_tissue.tsv"
    file = pd.read_csv(expression_report_loc, sep='\t', low_memory=False)
    df = pd.DataFrame(file)
    gene_names = df["Gene name"]

    if curr_gene == "Intergenic" or curr_gene == "intergenic":  # Expression level is unknown for Intergenic SNPs due to lack of mapping from GWAS to Human Protein Atlas
        return "-2\t"
    else:  # SNP is located in a gene
        if curr_gene in gene_names.values:
            express_range = df.loc[df["Gene name"] == curr_gene]

            if disease == "breast_carcinoma":
                express_range = express_range.loc[express_range["Tissue"] == "breast"]

                if "High" in express_range["Level"].values:
                    return expression_dict["High"] + "\t"
                elif "Medium" in express_range["Level"].values:
                    return expression_dict["Medium"] + "\t"
                elif "Low" in express_range["Level"].values:
                    return expression_dict["Low"] + "\t"
                elif "Not detected" in express_range["Level"].values:
                    return expression_dict["Not detected"] + "\t"
                else:
                    return "-1\t"

            if disease == "diabetes":
                express_range = express_range.loc[express_range["Tissue"] == "pancreas"]

                if "High" in express_range["Level"].values:
                    return expression_dict["High"] + "\t"
                elif "Medium" in express_range["Level"].values:
                    return expression_dict["Medium"] + "\t"
                elif "Low" in express_range["Level"].values:
                    return expression_dict["Low"] + "\t"
                elif "Not detected" in express_range["Level"].values:
                    return expression_dict["Not detected"] + "\t"
                else:
                    return "-1\t"

            if disease == "cardiovascular_disease":
                express_range = express_range.loc[express_range["Tissue"] == "heart muscle"]

                if "High" in express_range["Level"].values:
                    return expression_dict["High"] + "\t"
                elif "Medium" in express_range["Level"].values:
                    return expression_dict["Medium"] + "\t"
                elif "Low" in express_range["Level"].values:
                    return expression_dict["Low"] + "\t"
                elif "Not detected" in express_range["Level"].values:
                    return expression_dict["Not detected"] + "\t"
                else:
                    return "-1\t"

        else:  # Expression is not found at all
            return "-1\t"

def cleanup_single_gene(org_gene):
    cleaned_up_gene = org_gene.split(" ", 1)
    cleaned_up_gene = cleaned_up_gene[0]
    cleaned_up_gene = cleaned_up_gene.split("\\", 1)
    cleaned_up_gene = cleaned_up_gene[0]
    if (cleaned_up_gene[-1] == ",") or (cleaned_up_gene[-1] == ";"):
        cleaned_up_gene = cleaned_up_gene[:-1]
    return cleaned_up_gene

def gene_cleanup(genes):
    # gene clean-up, we will consider first mapped gene listed
    for i in range(len(genes)):
        org_gene = str(genes[i])
        genes[i] = cleanup_single_gene(org_gene)
    return genes
