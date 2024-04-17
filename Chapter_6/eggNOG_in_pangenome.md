# Let us match the eggNOG annotation to the pangenome results:
# Note this had to be done separately for E. coli and P. aeruginosa

Let's get the COGs and corresponding PROKKA IDs (to be able to merge with the Roary output):

```
# List all unique COGs:
 cut -f1,5 *.output.emapper.annotations | cut -f1 -d '@' | grep -v '#' > Prokka_id_COGs.txt
```

Now we need to match the PROKKA IDs to the Roary output to get the number of isolates - 
this will help identify which section of the pangenome the gene belongs to.
```
# first get a list of just PROKKA IDs:
cut -f1 Prokka_id_COGs.txt | sort | uniq > Prokka_ids.txt

# Now let's find them in the Roary output file: (gene_presence_absence.csv)
awk -F'"' 'NR==FNR{a[$1]; next} {for(i=1;i<=NF;i++) if($i in a) {print $i, $8}}' test_10_IDS.txt gene_presence_absence.csv > output.txt
```
Note we are using the '"' as a delimiter not ',' because some of the genes have commas in the names, but the values are separated by '"' as well, so we can ensure we extract the correct information.

Now we will merge the information to match the COGs and isolate numbers using a python script:
```
import pandas as pd
Query_gene = pd.read_csv("Key_gene_and_query.txt", names=["Query", "gene"], sep=' ')
query_isolate_no = pd.read_csv("Eggnog_prokka_id_isolate_number_NEW.csv", sep=' ', names=['Query', 'number'])
Query_cog = pd.read_csv("E_coli_query_cogs.txt", names=["Query", "cogs"], sep='\t')
Data = pd.merge(right=Query_gene, left=Query_cog, on='Query')
Data = Data.drop_duplicates()
cog_count = pd.merge(query_isolate_no, Data, on='Query')
cog_count = cog_count.drop('Query', axis=1)
cog_count = cog_count.drop('gene', axis=1)
cog_count = cog_count.drop_duplicates()
cog_count.to_csv("Cogs_by_isolate_count.csv", index=False)
```
To get a unique list of COGs and COG categories:
```
# list all unique COGs:
 cut -f5 *.output.emapper.annotations | cut -f1 -d '@' | grep -v '#' | grep -v 'eggNOG_OGs' > all_cogs.txt

# get a file of all corresponding COG categories:
cut -f7 *.output.emapper.annotations | grep -v '#' | grep -v 'COG_category' > all_cog_cat.txt

# after checking they are the same length, merge the files:

wc -l all_cog_cat.txt      # check number of lines in the file 
wc -l all_cogs.txt         # check number of lines in the file

# merge files:
paste -d ',' all_cogs.txt all_cog_cat.txt > Eggnog_cog_bycat.csv

# remove duplicates:
sort Eggnog_cog_bycat.csv | uniq > Eggnog_cog_bycat_uniq.csv
```

Some of these lines may have a COG corresponding to multiple COG categories, this will cause problems with our analysis (i.e. grouping by COG category)
Therefore, we will duplicate lines with multiple COG categories to have one COG category per line. 
For example:
COG  AK
goes to:
COG  A
COG  K

To achieve this we used the following Python code:
import pandas as pd

# eggnog COGs with cog cat
```
df = pd.read_csv("COG_cat_by_count_E_coli.csv")

# Function to split the letters in a string into a list
def split_letters(s):
    return list(s)

# Apply the function to the 'cog_cat' column and explode the DataFrame
df['cog_cat'] = df['cog_cat'].apply(split_letters)
df = df.explode('cog_cat', ignore_index=True)

# Display the modified DataFrame
print(df)

```

We now want to see which COG categories are present in the different pangenome sections:
```
# Continuing in the same python script
# Define the conditions for the new column
conditions = [
    (df['number'] > 9488.16),
    (df['number'] > 9104.8) & (df['number'] <= 9488.16),
    (df['number'] > 1437.6) & (df['number'] <= 9104.8)
]

# Define the corresponding values for each condition
values = ['core genome', 'softcore', 'shell']

# Add the new column based on the conditions
df['genome_type'] = np.select(conditions, values, default='cloud genome')

```

Now we can use these files to make the Figures






