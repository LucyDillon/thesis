# Lets match the RGI output to the pangenome results.
- we will use the PROKKA IDs to match the gene_presence_absence.csv (roary output) and we will use the RGI *.txt output files.
- First we need a unique list of RGI genes and PROKKA IDs, to do this I used the following bash commands (linux command line):

Let's get the AMR genes and corresponding PROKKA IDs (to be able to merge with the Roary output):

```
# List all unique genes and PROKKA id:
 cut -f1,9 *.txt | grep -v '#' > Prokka_id_genes.txt
```

Now we need to match the PROKKA IDs to the Roary output to get the number of isolates - 
this will help identify which section of the pangenome the gene belongs to.
```
# first get a list of just PROKKA IDs:
cut -f1 Prokka_id_genes.txt | sort | uniq > Prokka_ids.txt

# Now let's find them in the Roary output file: (gene_presence_absence.csv)
awk -F'"' 'NR==FNR{a[$1]; next} {for(i=1;i<=NF;i++) if($i in a) {print $i, $8}}' Prokka_ids.txt gene_presence_absence.csv > Query_isolate_no.txt
```
Note we are using the '"' as a delimiter not ',' because some of the genes have commas in the names, but the values are separated by '"' as well, so we can ensure we extract the correct information.

Now we will merge the information to match the genes and isolate numbers using a python script:
```
import pandas as pd
Query_gene = pd.read_csv("Key_gene_and_query.txt", names=["Query", "gene"], sep=' ')
query_isolate_no = pd.read_csv("Query_isolate_no.txt", sep=' ', names=['Query', 'number'])
Query_ARG = pd.read_csv("E_coli_query_ARGs.txt", names=["Query", ""], sep='\t')
Data = pd.merge(right=Query_gene, left=Query_ARG, on='Query')
Data = Data.drop_duplicates()
gene_count = pd.merge(query_isolate_no, Data, on='Query')
gene_count = gene_count.drop('Query', axis=1)
gene_count = gene_count.drop('gene', axis=1)
gene_count = gene_count.drop_duplicates()
gene_count.to_csv("genes_by_isolate_count.csv", index=False)
```
To get a unique list of genes and drug classes:
```
# list all unique RGI genes: where *txt is the RGI output file
 cut -f9 *.txt | grep -v 'Best_Hit_ARO'  > all_genes.txt

# get a file of all corresponding drugclasses:
cut -f15 *.txt | grep -v 'Drug Class' > all_drug_classes.txt

# after checking they are the same length, merge the files:

wc -l all_genes.txt      # check number of lines in the file 
wc -l all_drug_classes.txt         # check number of lines in the file

# merge files:
paste -d ',' all_genes.txt  all_drug_classes.txt > genes_drugclass.csv

# remove duplicates:
sort genes_drugclass.csv | uniq > genes_drugclass_uniq.csv
```

Some of these lines may have a gene corresponding to multiple drug classes, this will cause problems with our analysis (i.e. grouping by drug class)
Therefore, we will duplicate lines with multiple drug classes to have one drug class per line. 
For example:
gene1  aminoglycoside;carbapenem
goes to:
gene1  aminoglycoside
gene1  carbapenem

To achieve this we used the following Python code:
```
import pandas as pd
isolate_number = pd.read_csv("genes_drugclass_uniq.csv", names=['ORF_ID', 'isolate_no']) # where ORF id is the PROKKA id
RGI = pd.read_csv("RGI_drug_class_info_fixed.txt", sep='\t')
merged_data = pd.merge(RGI, isolate_number, on='ORF_ID', how='inner')
merged_data['isolate_no'].astype(str)
merged_data['Drug Class'].astype(str)
merged_data = merged_data.dropna()
split_data = []

for index, row in merged_data.iterrows():
    drug_classes = row['Drug Class'].split('; ')
    for drug_class in drug_classes:
        split_data.append([row['ORF_ID'], row['Best_Hit_ARO'], drug_class, row['isolate_no']])

split_df = pd.DataFrame(split_data, columns=merged_data.columns)

split_df_reduced = split_df[(split_df['Drug Class'] == 'cephalosporin')|(split_df['Drug Class'] == 'penam')| (split_df['Drug Class'] == 'fluoroquinolone antibiotic')| 
                            (split_df['Drug Class'] == 'tetracycline antibiotic')| (split_df['Drug Class'] == 'phenicol antibiotic')|(split_df['Drug Class'] == 'aminoglycoside antibiotic') 
                            | (split_df['Drug Class'] == 'macrolide antibiotic')| (split_df['Drug Class'] == 'monobactam')| (split_df['Drug Class'] == 'carbapenem')| (split_df['Drug Class'] == 'penem')|
                            (split_df['Drug Class'] == 'nitrofuran antibiotic')| (split_df['Drug Class'] == 'glycopeptide antibiotic')]

conditions = [
    (split_df_reduced['isolate_no'] > 9488.16),
    (split_df_reduced['isolate_no'] > 9104.8) & (split_df_reduced['isolate_no'] <= 9488.16),
    (split_df_reduced['isolate_no'] > 1437.6) & (split_df_reduced['isolate_no'] <= 9104.8)
]

# Define the corresponding values for each condition
values = ['core genome', 'softcore', 'shell']

# Add the new column based on the conditions
split_df_reduced['genome_type'] = np.select(conditions, values, default='cloud genome')
split_df_reduced.head()

```




