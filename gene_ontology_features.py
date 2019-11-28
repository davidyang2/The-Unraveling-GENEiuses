import pandas
import os

def add_to_header():
    return "GO FUNCTION\t"

def add_feature_for_curr_snp(disease, curr_gene):
    #print("Searching for " + str(curr_gene))
    cluster_report_loc = ""
    if disease == "diabetes":
        cluster_report_loc = "/Users/kavya/JHU/comp_bio/project/diabetes_functional_cluster.tsv"
    
    # We will now search the functional cluster tsv file to determine which function cluster the
    # gene belongs to. Note that many genes appear in multiple clusters, so we will pick the cluster
    # the gene is in that has the highest enrichment score.

    curr_gene = cleanup_single_gene(str(curr_gene))
    cluster_i = 1
    found_gene = False
    with open(cluster_report_loc, 'r') as file:
        read_a_new_cluster = file.readline().rstrip()
        while (read_a_new_cluster):
            #print("Looking at cluster " + str(cluster_i))
            #on_same_cluster = True
            table_name = read_a_new_cluster
            table_headers = file.readline().rstrip().split('\t')

            # keep looping through current cluster's table until we reach empty line
            on_same_cluster = file.readline().rstrip()
            while (on_same_cluster):
                data_row = on_same_cluster.split('\t')
                data_genes_row = data_row[5]
                data_genes_row = data_genes_row.rstrip().split(',')
                #print(data_genes_row)
                for d in data_genes_row:
                    d = d.rstrip().lstrip()
                    #print(d)
                    if str(d) == curr_gene:
                        found_gene = True
                        break

                if found_gene:
                    break
                on_same_cluster = file.readline().rstrip()

            if found_gene:
                break
            # if we read another newline / empty line, we have reached EOF
            read_a_new_cluster = file.readline().rstrip()
            #print(read_a_new_cluster)
            cluster_i = cluster_i + 1


    if (not found_gene):
        cluster_i = 0       # return 0 as label if we have not found the gene in any cluster

    #print(str(cluster_i))
    return str(cluster_i) + "\t"

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

def make_gene_list(disease, pos_loc, neg_loc, pos_indices, neg_indices):
    file = pandas.read_csv(pos_loc, sep='\t', lineterminator='\r')
    #pos_genes = file.MAPPED_GENE[pos_indices].values
    #pos_genes = pandas.read_csv(pos_loc, index_col = "MAPPED GENE", sep='\t', lineterminator='\r')
    df = pandas.DataFrame(file)
    pos_genes = df[df.columns[13]]
    pos_genes = pos_genes[pos_indices].values
    pos_genes = gene_cleanup(pos_genes)

    file = pandas.read_csv(neg_loc, sep='\t', lineterminator='\r', low_memory=False)
    df = pandas.DataFrame(file)
    neg_genes = df[df.columns[13]]
    neg_genes = neg_genes[neg_indices].values
    #neg_genes = file.MAPPED_GENE[neg_indices].values
    neg_genes = gene_cleanup(neg_genes)

    genes = set(pos_genes)
    genes.update(set(neg_genes))

    out_file_name = os.path.join(os.path.split(pos_loc)[0], disease + '_gene_list.txt')
    with open(out_file_name, "w+") as out_file:
        for g in genes:
            out_file.write(str(g) + "\n")