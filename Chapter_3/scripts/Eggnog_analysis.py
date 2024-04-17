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


Gene_families = eggnog_data['Gene_families']
Gene_families = Gene_families.sort_values(ascending=True)
Gene_families = Gene_families.drop_duplicates()
Gene_families = Gene_families.dropna()
Gene_families = Gene_families.astype(str)

antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
               'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erthromycin', 'fosfomycin',
               'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
               'tigecycline', 'tobramycin']

# Loop over the antibiotics
for antibiotic in antibiotics:
    # Make output files:
    antibiotic_table = final_table[final_table["antibiotic"] == antibiotic]
    antibiotic_table = antibiotic_table.drop_duplicates()
    antibiotic_table = antibiotic_table.set_index('genome_id')
    antibiotic_table = antibiotic_table.drop('antibiotic', axis=1)
    antibiotic_array = antibiotic_table.to_string(index=False, header=False)
    antibiotic_array_sep = antibiotic_array.replace(' ', ',')
    antibiotic_array_sep = antibiotic_array_sep.replace(',,,', ',')
    antibiotic_array_sep = antibiotic_array_sep.replace(',,', ',')

    # Set the outfile name based on the antibiotic
    outfile_name = antibiotic + '_Eggnog.arff'

    # Open the attribute output file
    att_out = open(outfile_name, "w")

    # Write the attribute list into the file
    i = 'GeneFamily'
    att_out.write('@RELATION    ' + i +'\n')
    for j in Gene_families:
        att_out.write('@ATTRIBUTE ' + j + '    REAL\n')
    name = 'phenotype'
    att_out.write('@ATTRIBUTE ' + name + '        {Susceptible, Resistant}\n')
    att_out.write('@DATA\n')
    # Write the array data to the file and close it
    att_out.write(antibiotic_array_sep)
    att_out.close()
