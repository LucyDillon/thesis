# Import packages:
import pandas as pd
import numpy as np

# Read in MIC & EUCAST data files:
Species_ids = pd.read_csv("All_patric_genus_species.csv")
EUCAST_MIC = pd.read_csv("EUCAST_MIC_values.csv")
MIC_data = pd.read_csv("ALL_patric_MIC_genomes.csv")
drug_class = pd.read_csv("merge_drugclass.csv")

# Read in all genomes:
eggnog_data = pd.read_csv("eggnog_data_all_genomes.csv")
eggnog_data.dropna()

eggnog_data['genome_id'] =eggnog_data['genome_id'].astype(str)

# Group data:
grouped = eggnog_data.groupby(["genome_id", "Gene_families"]).size().reset_index(name="Family_freq")

# Convert into a table:
table = pd.pivot_table(grouped, values=['Family_freq'], index=['genome_id'],
                    columns=['Gene_families'], fill_value=0)

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

phenotype['genome_id'] = phenotype['genome_id'].astype(str)

table = table.reset_index()

table.columns = table.columns.droplevel(0)
table.rename(columns = {list(table)[0]: 'genome_id'}, inplace = True)

#merge phenotype data with eggnog data
final_table = pd.merge(table, phenotype, on='genome_id')
# only keep rows with susceptible or res phenotype:
final_table = final_table[(final_table["add_phenotype"] == 'Susceptible') | (final_table["add_phenotype"] == 'Resistant')]

final_table['genome_id'] = final_table['genome_id'].astype(str)
Species_ids['genome_id'] = Species_ids['genome_id'].astype(str)

taxa_table = pd.merge(left=final_table, right=Species_ids, on='genome_id', how='left')

Gene_families = eggnog_data['Gene_families']
Gene_families = Gene_families.sort_values(ascending=True)
Gene_families = Gene_families.drop_duplicates()
Gene_families = Gene_families.dropna()
Gene_families = Gene_families.astype(str)


def process_ecoli_data(taxa_table, i, j):
    # Filter Escherichia coli strains and specified antibiotics
    E_coli = taxa_table[taxa_table['species'] == 'Escherichia coli']
    E_coli = E_coli[(E_coli['antibiotic'] == i)|(E_coli['antibiotic'] == j)]
    
    # Remove duplicates
    E_coli = E_coli.drop_duplicates()
    E_coli_reduced = E_coli[['genome_id', 'antibiotic', 'add_phenotype']]
    E_coli_reduced = E_coli_reduced.drop_duplicates()
    
    # Map 'add_phenotype' to binary values
    E_coli_reduced['phenotype'] = E_coli_reduced['add_phenotype'].map({'Susceptible': 0, 'Resistant': 1})
    
    # Pivot table for antibiotic-phenotype relationship
    E_coli_AB = pd.pivot_table(E_coli_reduced, values=['phenotype'], index=['genome_id'], columns=['antibiotic'])
    E_coli_AB = E_coli_AB.reset_index()
    E_coli_AB.columns = E_coli_AB.columns.droplevel(0)
    E_coli_AB.rename(columns={list(E_coli_AB)[0]: 'genome_id'}, inplace=True)
    
    # Separate data by phenotypes
    SR = E_coli_AB[(E_coli_AB[i[0]] == 0) & (E_coli_AB[j[1]] == 1)]
    SR.insert(0, 'phenotype', 'SR')
    RS = E_coli_AB[(E_coli_AB[i[0]] == 1) & (E_coli_AB[j[1]] == 0)]
    RS.insert(0, 'phenotype', 'RS')
    RR = E_coli_AB[(E_coli_AB[i[0]] == 1) & (E_coli_AB[j[1]] == 1)]
    RR.insert(0, 'phenotype', 'RR')
    SS = E_coli_AB[(E_coli_AB[i[0]] == 0) & (E_coli_AB[j[1]] == 0)]
    SS.insert(0, 'phenotype', 'SS')
    
    # Concatenate all phenotypes
    E_coli_by_pheno = pd.concat([SS, SR, RS, RR])
    E_coli_by_pheno_reduced = E_coli_by_pheno.drop(i, axis=1)
    E_coli_by_pheno_reduced = E_coli_by_pheno_reduced.drop(j, axis=1)
    # Prepare gene hits data
    E_coli_genehits = E_coli.drop(['antibiotic', 'add_phenotype', 'species', 'EUCAST_species', 'Genus'], axis=1)
    del E_coli_genehits[E_coli_genehits.columns[-1]]
    
    # Merge data to get the final DataFrame
    New_E_coli = pd.merge(left=E_coli_genehits, right=E_coli_by_pheno_reduced, on='genome_id', how='right')
    
    # Drop duplicates and set index
    New_E_coli = New_E_coli.drop_duplicates()
    New_E_coli = New_E_coli.set_index('genome_id')
    
    # Convert DataFrame to ARFF format
    E_coli_array = New_E_coli.to_string(index=False, header=False)
    E_coli_array_sep = E_coli_array.replace(' ', ',')
    E_coli_array_sep = E_coli_array_sep.replace(',,,', ',')
    E_coli_array_sep = E_coli_array_sep.replace(',,', ',')
        
    att_out = open(i + '_' + j + '_MDR_TEST_TEST.arff', "w")

    # Write the attribute list into the file
    x = 'GeneFamily'
    att_out.write('@RELATION    ' + x +'\n')
    for y in Gene_families:
        att_out.write('@ATTRIBUTE ' + y + '    REAL\n')
    name = 'phenotype'
    att_out.write('@ATTRIBUTE ' + name + '        {SS, SR, RS, RR}\n')
    att_out.write('@DATA\n')
    # Write the array data to the file and close it
    att_out.write(E_coli_array_sep)
    att_out.close() 
    
# Usage example
Antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
               'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erthromycin', 'fosfomycin',
               'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
               'tigecycline', 'tobramycin']

for i in Antibiotics:
    for j in Antibiotics:
        process_ecoli_data(taxa_table, i, j)
        