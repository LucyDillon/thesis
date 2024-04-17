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

# MIC analysis:
MIC_reduced = MIC_data[(MIC_data["phenotype"] == "Resistant") | (MIC_data["phenotype"]== "Susceptible")]
data_reduced = data.drop('Gene', axis =1)
data_antibioticed = data_reduced.antibioticby(['genome_id', 'RGIDrugClass'])['GeneHit'].sum().reset_index()
# If gene hit > 0 genome is resistant, if 1 < susceptible
data_antibioticed['pheno'] = np.where(data_antibioticed['GeneHit'] > 0, "Resistant", "Susceptible")
# Now merge with the drug class to get the matching AB to the phenotype
RGI_data = pd.merge(data_antibioticed, drug_class, on='RGIDrugClass')
# Reduce columns:
RGI_data_reduced = RGI_data[['genome_id', 'GeneHit', 'pheno', 'antibiotic']]
RGI_data_reduced = RGI_data_reduced.dropna()

#Merge Species_ids and EUCAST_MIC with common column "EUCAST_species"
EUCAST_species_MIC = pd.merge(EUCAST_MIC, Species_ids, on="EUCAST_species")

#Merge EUCAST_species_MIC with original MIC values 
Overall_EUCAST_data = pd.merge(left=EUCAST_species_MIC, right=MIC_reduced, how='left', on=['genome_id', 'antibiotic'])
Overall_EUCAST_data = Overall_EUCAST_data.dropna(subset=['phenotype'])

#Take average of MIC values by antibioticing by genome id and antibiotic
average_Overall_EUCAST_data = Overall_EUCAST_data.antibioticby(['genome_id', 'antibiotic']).agg({'MIC': ['mean']})
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
# Now need to merge with the phenotype 
all_data = pd.merge(RGI_data_reduced, phenotype, on=['genome_id','antibiotic'])
# only keep rows with susceptible or res phenotype:
final_table = all_data[(all_data["add_phenotype"] == 'Susceptible') | (all_data["add_phenotype"] == 'Resistant')]
antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
               'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erthromycin', 'fosfomycin',
               'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
               'tigecycline', 'tobramycin']

for antibiotic in antibiotics:
    # Filter the DataFrame for the specific antibiotic
    antibiotic_df = final_table[final_table['antibiotic'] == antibiotic]
    antibiotic_df = antibiotic_df.drop_duplicates()
    
    # Filter rows for predicted 'Susceptible' and 'Resistant'
    predicted_sus = antibiotic_df[antibiotic_df['pheno'] == 'Susceptible']
    predicted_res = antibiotic_df[antibiotic_df['pheno'] == 'Resistant']
    
    # Calculate true positives, true negatives, false positives, and false negatives
    true_positive = predicted_sus[predicted_sus['add_phenotype'] == 'Susceptible'].shape[0]
    true_negative = predicted_res[predicted_res['add_phenotype'] == 'Resistant'].shape[0]
    false_positive = predicted_res[predicted_res['add_phenotype'] == 'Susceptible'].shape[0]
    false_negative = predicted_sus[predicted_sus['add_phenotype'] == 'Resistant'].shape[0]
    
    # Create the confusion matrix
    confusion_matrix = pd.DataFrame(
        [[true_negative, false_positive], [false_negative, true_positive]],
        index=['Resistant', 'Susceptible'],
        columns=['Resistant', 'Susceptible']
    )

    print(f"Confusion matrix for {antibiotic}:\n{confusion_matrix}\n")
    
