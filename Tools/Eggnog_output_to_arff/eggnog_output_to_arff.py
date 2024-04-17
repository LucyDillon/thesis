# Import packages:
import pandas as pd
import numpy as np
import argparse

# Add menu for command line usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run eggnog to arff tool parameters')
    parser.add_argument('-file', action='store', dest='file', required=True,
                        help='eggnog output file')
    parser.add_argument('-o', action='store', dest='output_file',
                        help='Name of output file ext .arff', required=True)
    parser.add_argument('-c_attr_name', action='store', dest='cat_attr_name',
                        help='Catagorical variable column title', required=True)
    parser.add_argument('-catagory1', action='store', dest='catagory1',
                        help='Catagory 1 e.g. Susceptible', required=True)
    parser.add_argument('-catagory2', action='store', dest='catagory2',
                        help='Catagory 2 e.g. Resistant', required=True)
    parser.add_argument('-catagory3', action='store', dest='catagory3',
                        help='Catagory 3 e.g. Intermediate', required=False)

options = parser.parse_args()
#Read in eggnog file:
eggnog = pd.read_csv(options.file, sep='\\t', comment='#',names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'] )
# Subset table only selecting column 1 and 5
eggnog_cogs = eggnog[['query', 'eggNOG_OGs']]

eggnog_cogs['query'] = eggnog_cogs['query'].str.split('_').str[0]
eggnog_cogs['eggNOG_OGs'] = eggnog_cogs['eggNOG_OGs'].apply(lambda s: s.split('@')[0])

eggnog_cogs = eggnog_cogs.rename(columns={'query': 'genome_id', 'eggNOG_OGs': 'cogs'})

grouped = eggnog_cogs.groupby(["genome_id", "cogs"]).size().reset_index(name="cog_freq")

table = pd.pivot_table(grouped, values=['cog_freq'], index=['genome_id'],
                    columns=['cogs'], fill_value=0)

#add empty phenotype column:
# this is now catagorical vairable:
table[options.cat_attr_name] = np.nan

Cogs =eggnog_cogs['cogs']
Cogs= Cogs.sort_values(ascending=True)
Cogs = Cogs.drop_duplicates()

table_array = table.to_string(index=False, header=False)
table_array_sep = table_array.replace(' ', ',')
table_array_sep = table_array_sep.replace(',,,', ',')
table_array_sep = table_array_sep.replace(',,', ',')
outfile = (options.output_file)
att_out = open(outfile, "w")
# Write an attribute list into the first part of the file:
i = ('Cog')
att_out.write('@RELATION    ' + i +'\n')
for j in Cogs:
    att_out.write("@ATTRIBUTE " + i + "-" + j + "    REAL\n")
name = (options.cat_attr_name)
att_out.write('@ATTRIBUTE ' + name + '		{' + options.catagory1 + ', ' + options.catagory2 + '}\n')
att_out.write('@DATA\n')
# Write the array data to the bottom of the file and close the file:
att_out.write(table_array_sep)
att_out.write('\n')
att_out.close()
