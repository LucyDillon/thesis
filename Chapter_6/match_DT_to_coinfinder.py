# This file is an example of how we matched the Ciprofloxacin and Ceftriaxone

import pandas as pd

# Read the gene and coinfinder networks
genes_df = pd.read_csv("Cef_cip_genes.txt", header=None, names=["Gene"])
csv_df = pd.read_csv("coinfinder_assoc_cogs.csv") # also used this script for disassociated networks. I.e. switch the csv file 

# Initialize an empty list to store gene pairs
gene_pairs_list = []

# Iterate through each gene in the source column
for gene in genes_df["Gene"]:
    # Filter rows where the gene is present in the source column
    source_rows = csv_df[csv_df["COG_source"].str.contains(gene)]

    # Check if any of the other genes are in the target column
    for target_gene in genes_df["Gene"]:
        gene_pair_rows = source_rows[source_rows["COG_target"].str.contains(target_gene)]
        if not gene_pair_rows.empty:
            gene_pairs_list.append(gene_pair_rows[["COG_source", "COG_target"]])

# Concatenate the reduced DataFrames into a single DataFrame
result_df = pd.concat(gene_pairs_list, ignore_index=True)

# Save the result DataFrame to a CSV file
result_df.to_csv("gene_pairs_associated_Cef_cip.csv", index=False)
