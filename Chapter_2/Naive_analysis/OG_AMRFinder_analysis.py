# Import packages:
import pandas as pd
import numpy as np
# Read in MIC & EUCAST data files:
Species_ids = pd.read_csv("All_patric_genus_species.csv")
EUCAST_MIC = pd.read_csv("EUCAST_MIC_values.csv")
MIC_data = pd.read_csv("ALL_patric_MIC_genomes.csv")
drug_class = pd.read_csv("merge_drugclass.csv")
AMRFinder_drugclass = pd.read_csv("AMR_genes_DrugClass.csv", names=['Gene', 'AMRFinderDrugClass'])
# Read in all genomes:
complete_genomes = pd.read_csv("AMRFinder_R_format_complete.csv")
AMRFinder_WGS = pd.read_csv("myRinputfile.csv", names=['genome_id', 'Gene', 'GeneHit'])
complete_genomes = complete_genomes[['genome_id', 'Gene', 'GeneHit', 'AMRFinderDrugClass']]
AMRFinder_WGS_drugclass = pd.merge(left=AMRFinder_WGS, right=AMRFinder_drugclass, on='Gene')
Species = Species_ids[['genome_id', 'species']]
taxa_info = MIC_data[['genome_id', 'Phylum', 'Class', 'Order', 'Family']]

# Merge resfinder WGS and complete datasets:
AMRFinder = [AMRFinder_WGS_drugclass, complete_genomes]
AMRFinder_all = pd.concat(AMRFinder)
# MIC analysis:
MIC_reduced = MIC_data[(MIC_data["phenotype"] == "Resistant") | (MIC_data["phenotype"]== "Susceptible")]


#Merge Species_ids and EUCAST_MIC with common column "EUCAST_species"
EUCAST_species_MIC = pd.merge(EUCAST_MIC, Species_ids, on="EUCAST_species")

#Merge EUCAST_species_MIC with original MIC values 
Overall_EUCAST_data = pd.merge(left=EUCAST_species_MIC, right=MIC_reduced, how='left', on=['genome_id', 'antibiotic'])
Overall_EUCAST_data = Overall_EUCAST_data.dropna(subset=['phenotype'])

#Take average of MIC values by grouping by genome id and antibiotic
average_Overall_EUCAST_data = Overall_EUCAST_data.groupby(['genome_id', 'antibiotic']).agg({'MIC': ['mean']})
average_Overall_EUCAST_data.columns = ['MIC_mean']
average_Overall_EUCAST_data = average_Overall_EUCAST_data.reset_index()

#merge the average MIC with the overall data to get the extra columns
test_data = pd.merge(average_Overall_EUCAST_data, Overall_EUCAST_data, how='left', left_on=['genome_id', 'antibiotic'], right_on=['genome_id', 'antibiotic'])

#Make a reduced overall dataframe called "Reduced_EUCAST_df"
Reduced_EUCAST_df = test_data[["genome_id", "antibiotic", "MIC_mean", "EUCAST_MIC", "EUCAST_species"]]

#Separate species names and phenotype by separating by underscore delimeter into 2 separate columns
Reduced_EUCAST_df[['EUCAST_species','phenotype']] = Reduced_EUCAST_df['EUCAST_species'].str.split('_',expand=True)

#Make two subset tables resistant_EUCAST and susceptible_EUCAST to manipulate
resistant_EUCAST = Reduced_EUCAST_df[Reduced_EUCAST_df["phenotype"] == "Res"]
susceptible_EUCAST = Reduced_EUCAST_df[Reduced_EUCAST_df["phenotype"] == "sus"]

Overall_data = pd.merge(resistant_EUCAST, susceptible_EUCAST, how='left', on=['genome_id', 'antibiotic', 'MIC_mean', 'EUCAST_species'])

def add_phenotype(MIC_mean, Res, Sus):
        if MIC_mean > Res: 
                return("Resistant")
        elif MIC_mean <= Sus:
                return("Susceptible")
        else: 
                return("Intermediate")

Overall_data['add_phenotype'] =Overall_data.apply (lambda row : add_phenotype(row['MIC_mean'], row['EUCAST_MIC_x'], row['EUCAST_MIC_y']), axis = 1)

#merge overall data and drug class table (this is so it can be merged with RGI on genome id and antibiotic)
all_EUCAST_data = pd.merge(Overall_data, drug_class, on="antibiotic")
# Reduce all_EUCAST_data:
all_EUCAST_data_final = all_EUCAST_data[['genome_id', 'antibiotic', 'add_phenotype']]
phenotype = all_EUCAST_data_final[["antibiotic", "genome_id", "add_phenotype"]]


#merge phenotype data with RGI data
data_table = pd.merge(AMRFinder_all, drug_class, on=['AMRFinderDrugClass'])

final_table = pd.merge(data_table, phenotype, on=["genome_id", "antibiotic"])

final_table_reduced = final_table[["genome_id", "Gene", "GeneHit", "antibiotic", "add_phenotype"]]

# only keep rows with susceptible or res phenotype:
final_table = final_table_reduced[(final_table_reduced["add_phenotype"] == 'Susceptible') | (final_table_reduced["add_phenotype"] == 'Resistant')]

Pivot_table = pd.pivot_table(final_table, values='GeneHit', index='genome_id', columns='Gene', fill_value=0)
Pivot_table['sum'] = Pivot_table[list(Pivot_table.columns)].sum(axis=1)
Pivot_table['pheno'] = np.where(Pivot_table['sum'] > 0, "Resistant", "Susceptible")
Pivot_table = Pivot_table[["pheno"]]

phenotype =phenotype.set_index('genome_id')
#merge phenotype data with eggnog data
final_table = Pivot_table.join(phenotype)

final_table = final_table[(final_table["add_phenotype"] == 'Susceptible') | (final_table["add_phenotype"] == 'Resistant')]

# Split the dataframe into groups by 'antibiotic'
groups = final_table.groupby('antibiotic')

# Create an empty dictionary to store the confusion matrices for each group
confusion_matrices = {}

# Iterate over each group and create a confusion matrix
for group_name, group in groups:
    # Create a boolean mask for true positives and true negatives
    tp_mask = (group['pheno'] == 'Susceptible') & (group['add_phenotype'] == 'Susceptible')
    tn_mask = (group['pheno'] == 'Resistant') & (group['add_phenotype'] == 'Resistant')
    fp_mask = (group['pheno'] == 'Susceptible') & (group['add_phenotype'] == 'Resistant')
    fn_mask = (group['pheno'] == 'Resistant') & (group['add_phenotype'] == 'Susceptible')
    
    # Calculate the number of true positives, true negatives, false positives, and false negatives
    tp = tp_mask.sum()
    tn = tn_mask.sum()
    fp = fp_mask.sum()
    fn = fn_mask.sum()
    
    # Create a confusion matrix for the group
    matrix = pd.DataFrame({
        'True Positive': [tp],
        'True Negative': [tn],
        'False Positive': [fp],
        'False Negative': [fn]
    }, index=[group_name])
    
    # Store the confusion matrix in the dictionary
    confusion_matrices[group_name] = matrix
    
# Concatenate the confusion matrices for all groups into a single dataframe
result = pd.concat(confusion_matrices.values(), axis=0)

# Print the result
print(result)
