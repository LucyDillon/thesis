# Import packages:
import pandas as pd
import numpy as np

# Read in MIC & EUCAST data files:
Species_ids = pd.read_csv("All_patric_genus_species.csv")
EUCAST_MIC = pd.read_csv("EUCAST_MIC_values.csv")
MIC_data = pd.read_csv("ALL_patric_MIC_genomes.csv")
drug_class = pd.read_csv("merge_drugclass.csv")
# Read in all genomes:
data = pd.read_csv("ALL_RGI_genomes.csv")

Species = Species_ids[['genome_id', 'Genus']]

taxa_info = MIC_data[['genome_id', 'Phylum', 'Class', 'Order', 'Family']]

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
data_table = pd.merge(data, drug_class, on=['RGIDrugClass'])

final_table = pd.merge(data_table, phenotype, on=["genome_id", "antibiotic"])

final_table_reduced = final_table[["genome_id", "Gene", "GeneHit", "antibiotic", "add_phenotype"]]

# only keep rows with susceptible or res phenotype:
final_table = final_table_reduced[(final_table_reduced["add_phenotype"] == 'Susceptible') | (final_table_reduced["add_phenotype"] == 'Resistant')]

# Convert into a table:
table = pd.pivot_table(final_table, values=['GeneHit'], index=['genome_id'],
                    columns=['Gene'], fill_value=0)

table = table.reset_index()
table.columns = table.columns.droplevel(0)
table.rename(columns = {list(table)[0]: 'genome_id'}, inplace = True)

# Add genus to the table:
taxonomy_table = pd.merge(left=table, right=Species, on='genome_id')
# drop duplicates:
taxonomy_table = taxonomy_info_2_1.drop_duplicates()

# Get phenotype info ready to merge with the table above:
pheno = final_table[['genome_id', 'antibiotic', 'add_phenotype']] 
pheno = pheno.drop_duplicates()

#merge phenotype data with RGI data
taxonomy_pheno_table= pd.merge(left=taxonomy_table, right=pheno, on='genome_id')

#Get a list of the unique genes:
Genes = taxonomy_info_2['Gene']
Genes = Genes.sort_values(ascending=True)
Genes = Genes.drop_duplicates()


#Make a arff file for each genus:            
def write_to_arff_genus(antibiotic, genu):
    antibiotic_df = taxonomy_pheno_table[taxonomy_pheno_table['antibiotic'] == antibiotic]
    antibiotic_df = taxonomy_pheno_table[taxonomy_pheno_table['Genus'] == genus]
    antibiotic_df = antibiotic_df.drop_duplicates()
    antibiotic_df = antibiotic_df.set_index('genome_id')
    antibiotic_df = antibiotic_df.drop('antibiotic', axis=1)
    antibiotic_df = antibiotic_df.drop('Genus', axis=1)
    antibiotic_df_array = antibiotic_df.to_string(index=False, header=False)
    antibiotic_df_array_sep = antibiotic_df_array.replace(' ', ',')
    antibiotic_df_array_sep = antibiotic_df_array_sep.replace(',,,', ',')
    antibiotic_df_array_sep = antibiotic_df_array_sep.replace(',,', ',')
    outfile_genus = f"{antibiotic}_{genus}_taxa_files.arff"
    att_out_genus = open(outfile_genus, "w")
    i = 'Gene'
    att_out_genus.write('@RELATION    ' + i + '\n')
    for x in Genes:
        att_out_genus.write('@ATTRIBUTE Gene-' + x + '  REAL\n')
    att_out_genus.write('@ATTRIBUTE phenotype		{Susceptible, Resistant}\n')
    att_out_genus.write('@DATA\n')
    att_out_genus.write(antibiotic_df_array_sep)
    att_out_genus.close()
                

antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
               'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erthromycin', 'fosfomycin',
               'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
               'tigecycline', 'tobramycin'] 

genuses = ['Achromobacter', 'Acinetobacter', 'Campylobacter', 'Citrobacter', 
           'Corynebacterium', 'Enterobacter', 'Escherichia', 'Klebsiella', 
           'Kluyvera', 'Listeria', 'Morganella', 'Neisseria', 'Proteus', 'Providencia', 
           'Raoultella', 'Salmonella', 'Serratia', 'Shigella', 'Streptococcus']

for antibiotic in antibiotics:
    for genus in genuses:
        write_to_arff_genus(antibiotic, genus)
        

# Make a arff file that is all but one genus:

def write_to_arff_except_genus(antibiotic, genus):
    antibiotic_df = taxonomy_pheno_table[taxonomy_pheno_table['antibiotic'] == antibiotic]
    antibiotic_df = taxonomy_pheno_table[taxonomy_pheno_table['Genus'] != genus]
    antibiotic_df = antibiotic_df.drop_duplicates()
    antibiotic_df = antibiotic_df.set_index('genome_id')
    antibiotic_df = antibiotic_df.drop('antibiotic', axis=1)
    antibiotic_df = antibiotic_df.drop('Genus', axis=1)
    antibiotic_df_array = antibiotic_df.to_string(index=False, header=False)
    antibiotic_df_array_sep = antibiotic_df_array.replace(' ', ',')
    antibiotic_df_array_sep = antibiotic_df_array_sep.replace(',,,', ',')
    antibiotic_df_array_sep = antibiotic_df_array_sep.replace(',,', ',')
    outfile_genus = f"{antibiotic}_{genus}_not_taxa_files.arff"
    att_out_genus = open(outfile_genus, "w")
    i = 'Gene'
    att_out_genus.write('@RELATION    ' + i + '\n')
    for x in Genes:
        att_out_genus.write('@ATTRIBUTE Gene-' + x + '  REAL\n')
    att_out_genus.write('@ATTRIBUTE phenotype		{Susceptible, Resistant}\n')
    att_out_genus.write('@DATA\n')
    att_out_genus.write(antibiotic_df_array_sep)
    att_out_genus.close()
                

antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
               'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erthromycin', 'fosfomycin',
               'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
               'tigecycline', 'tobramycin'] 

excluded_genuses = ['Achromobacter', 'Acinetobacter', 'Campylobacter', 'Citrobacter', 
                    'Corynebacterium', 'Enterobacter', 'Escherichia', 'Klebsiella', 
                    'Kluyvera', 'Listeria', 'Morganella', 'Neisseria', 'Proteus', 'Providencia', 
                    'Raoultella', 'Salmonella', 'Serratia', 'Shigella', 'Streptococcus']

for antibiotic in antibiotics:
    for genus in excluded_genuses:
        write_to_arff_except_genus(antibiotic, genus)
