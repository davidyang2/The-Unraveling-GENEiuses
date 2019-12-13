import pandas as pd
import os
import sys


def add_to_header(disease):
    if disease == "breast_carcinoma":
        return "BREAST_A_LEVEL\tBREAST_G_LEVEL\tBREAST_M_LEVEL\t"

    if disease == "cardiovascular_disease":
        return "HEART_EXPRESSION_LEVEL\t"

    if disease == "diabetes":
        return "PANCREAS_EG_LEVEL\tPANCREAS_IOL_LEVEL\tKIDNEY_CG_LEVEL\tKIDNEY_T_LEVEL\t"

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
        if disease == "diabetes":
            return "-2\t-2\t-2\t-2\t"
        if disease == "breast_carcinoma":
            return "-2\t-2\t-2\t"
        else:
            return "-2\t"
    else:  # SNP is located in a gene
        if curr_gene in gene_names.values:
            express_range = df.loc[df["Gene name"] == curr_gene]

            if disease == "breast_carcinoma":
                # express_range = express_range.loc[express_range["Tissue"] == "breast"]

                print(curr_gene)
                expression_result = ""
                genes = list(df["Gene name"].values)
                tissues = list(df["Tissue"].values)
                level = list(df["Level"].values)
                for i in range(len(tissues)):
                    if genes[i] == curr_gene and tissues[i] == "breast":
                        expression_result += expression_dict[level[i]] + "\t"

                if expression_result == "":
                    return "-1\t-1\t-1\t"
                return expression_result

            if disease == "diabetes":
                # express_range_pancreas = express_range.loc[express_range["Tissue"] == "pancreas"]
                # express_range_kidney = express_range.loc[express_range["Tissue"] == "kidney"]

                expression_result = ""
                # print(curr_gene)

                genes = list(df["Gene name"].values)
                tissues = list(df["Tissue"].values)
                level = list(df["Level"].values)
                for i in range(len(tissues)):
                    if genes[i] == curr_gene and tissues[i] == "pancreas":
                        expression_result += expression_dict[level[i]] + "\t"
                for i in range(len(tissues)):
                    if genes[i] == curr_gene and tissues[i] == "kidney":
                        expression_result += expression_dict[level[i]] + "\t"

                if expression_result == "":
                    return "-1\t-1\t-1\t-1\t"
                return expression_result

            if disease == "cardiovascular_disease":
                # express_range = express_range.loc[express_range["Tissue"] == "heart muscle"]
                expression_result = ""
                genes = list(df["Gene name"].values)
                tissues = list(df["Tissue"].values)
                level = list(df["Level"].values)
                for i in range(len(tissues)):
                    if genes[i] == curr_gene and tissues[i] == "heart muscle":
                        expression_result += expression_dict[level[i]] + "\t"

                if expression_result == "":
                    return "-1\t"
                return expression_result

        else:  # Expression is not found at all
            if disease == "diabetes":
                return "-1\t-1\t-1\t-1\t"
            if disease == "breast_carcinoma":
                return "-1\t-1\t-1\t"
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
