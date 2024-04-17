import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix

MIC_data = pd.read_csv("ALL_patric_MIC_genomes.csv")
EUCAST_MIC = pd.read_csv("EUCAST_MIC_values.csv")
Species_ids = pd.read_csv("All_patric_genus_species.csv")
DrugClass_info = pd.read_csv("merge_drugclass.csv")

RGI_complete = pd.read_csv("RGI_R_format.csv")
RGI_WGS = pd.read_csv("RGI_input_final.csv", names=['genome_id', 'Gene', 'GeneHit'])
RGI_drugclass = pd.read_csv("RGI_genes_drugclass.csv", names=['Gene', 'RGIDrugClass'])
RGI_WGS_drugclass = pd.merge(left=RGI_WGS, right=RGI_drugclass, on='Gene', how='inner')
RGI_complete = RGI_complete[['genome_id', 'Gene', 'GeneHit', 'RGIDrugClass']]
RGI = [RGI_complete, RGI_WGS_drugclass]
RGI_all = pd.concat(RGI)

# Read in AMR Storfs:
RGI_storf = pd.read_csv("RGI_storf_input.csv", names=['genome_id', 'Gene', 'GeneHit'])
Storf_drugclass = pd.read_csv("unique_gene_drugclass.txt", sep='\t', names=['Gene', 'RGIDrugClass'])

RGI_storf_drugclass = pd.merge(RGI_storf, Storf_drugclass, on='Gene')
RGI_storf_all = [RGI_all,RGI_storf_drugclass]
RGI_ALL_AND_STORF = pd.concat(RGI_storf_all)

MIC_reduced = MIC_data[(MIC_data["phenotype"] == "Resistant") | (MIC_data["phenotype"] == "Susceptible")]

# Merge Species_ids and EUCAST_MIC with common column "EUCAST_species"
EUCAST_species_MIC = pd.merge(EUCAST_MIC, Species_ids, on="EUCAST_species")

# Merge EUCAST_species_MIC with original MIC values
Overall_EUCAST_data = pd.merge(left=EUCAST_species_MIC, right=MIC_reduced, how='left', on=['genome_id', 'antibiotic'])

Overall_EUCAST_data = Overall_EUCAST_data.dropna(subset=['phenotype'])

# Take average of MIC values by grouping by genome id and antibiotic
average_Overall_EUCAST_data = Overall_EUCAST_data.groupby(['genome_id', 'antibiotic']).agg({'MIC': ['mean']})
average_Overall_EUCAST_data.columns = ['MIC_mean']
average_Overall_EUCAST_data = average_Overall_EUCAST_data.reset_index()
# merge the average MIC with the overall data to get the extra columns
test_data = pd.merge(average_Overall_EUCAST_data, Overall_EUCAST_data, how='left', left_on=['genome_id', 'antibiotic'],
                     right_on=['genome_id', 'antibiotic'])
# Make a reduced overall dataframe called "Reduced_EUCAST_df"
Reduced_EUCAST_df = test_data[["genome_id", "antibiotic", "MIC_mean", "EUCAST_MIC", "EUCAST_species"]]
# Separate species names and phenotype by separating by underscore delimeter into 2 separate columns
Reduced_EUCAST_df[['EUCAST_species', 'phenotype']] = Reduced_EUCAST_df['EUCAST_species'].str.split('_', expand=True)
# Make two subset tables resistant_EUCAST and susceptible_EUCAST to manipulate
resistant_EUCAST = Reduced_EUCAST_df[Reduced_EUCAST_df["phenotype"] == "Res"]
susceptible_EUCAST = Reduced_EUCAST_df[Reduced_EUCAST_df["phenotype"] == "sus"]

Overall_data = pd.merge(resistant_EUCAST, susceptible_EUCAST, how='left', on=['genome_id', 'antibiotic',
                                                                              'MIC_mean', 'EUCAST_species'])

def add_phenotype(MIC_mean, Res, Sus):
    if MIC_mean > Res:
        return "Resistant"
    elif MIC_mean <= Sus:
        return "Susceptible"
    else:
        return "Intermediate"
    
    
Overall_data['add_phenotype'] = Overall_data.apply(lambda row: add_phenotype(row['MIC_mean'], row['EUCAST_MIC_x'],
                                                                             row['EUCAST_MIC_y']), axis=1)

# Read in amr drug class table to merge with overall data
drug_class = pd.read_csv("merge_drugclass.csv")
# merge overall data and drug class table (this is so it can be merged with RGI on genome id and antibiotic)
all_EUCAST_data = pd.merge(Overall_data, drug_class, on="antibiotic")

# Reduce all_EUCAST_data:
all_EUCAST_data_final = all_EUCAST_data[['genome_id', 'antibiotic', 'add_phenotype']]

# merge with EUCAST info to get normal AB
RGI_EUCAST = pd.merge(left=RGI_ALL_AND_STORF, right=DrugClass_info, on='RGIDrugClass')

RGI_EUCAST = RGI_EUCAST.drop_duplicates()

all_EUCAST_data_final['genome_id'] = all_EUCAST_data_final['genome_id'].astype(str)
RGI_EUCAST['genome_id'] = RGI_EUCAST['genome_id'].astype(str)

# Merge RGI data and EUCAST data:
FINAL_DATA = pd.merge(left=all_EUCAST_data_final, right=RGI_EUCAST, on=['genome_id', 'antibiotic'])

# reduce to essential columns only and drop duplicate data:
RGI_final = FINAL_DATA[['genome_id', 'antibiotic', 'add_phenotype', 'Gene', 'GeneHit']]


RGI_final = RGI_final[(RGI_final["add_phenotype"] == 'Susceptible') | (RGI_final["add_phenotype"] == 'Resistant')]



def antibiotic_classification(antibiotic):
    antibiotic = RGI_final[RGI_final['antibiotic'] == antibiotic]
    antibiotic = antibiotic[['genome_id', 'add_phenotype', 'GeneHit']]
    antibiotic = antibiotic[(antibiotic["add_phenotype"] == "Resistant") |
                            (antibiotic["add_phenotype"] == "Susceptible")]
    antibiotic_grouped = antibiotic.groupby(['genome_id', 'add_phenotype']).size().reset_index(name="all_genes")
    x = antibiotic_grouped['all_genes'].values.reshape(-1, 1)
    y = antibiotic_grouped['add_phenotype'].values
    xtrain, xtest, ytrain, ytest = train_test_split(
        x, y, test_size=0.25, random_state=0)
    sc_x = StandardScaler()
    xtrain = sc_x.fit_transform(xtrain)
    xtest = sc_x.transform(xtest)
    classifier = LogisticRegression(random_state=0)
    classifier.fit(xtrain, ytrain)
    print(classifier.coef_)
    print(classifier.intercept_)
    y_pred = classifier.predict(xtest)
    cm = confusion_matrix(ytest, y_pred)
    print("Confusion Matrix : \n", cm)
    print("Accuracy : ", accuracy_score(ytest, y_pred))

    
ab_list = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'cefepime', 'ceftriaxone', 'chloramphenicol',
           'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erythromycin', 'fosfomycin',
           'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline',
           'tigecycline', 'tobramycin']


for ab in ab_list:
    print(ab)
    antibiotic_classification(ab)
