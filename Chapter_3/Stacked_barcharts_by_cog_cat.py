# Import packages:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read in the data:
data = pd.read_csv("Gene_families_presence_absence.csv")

# Remove absent data:
data_reduced = data[data["presence_absence"] == 1]
data = data[data["Antibiotic"] != 'moxifloxacin']
data = data[data["Antibiotic"] != 'fosfomycin']
data = data[data["Antibiotic"] != 'amoxicillin']

# Stacked barchart of the 3 main categories:
data_grouped_overall = data_reduced.groupby(["Overall_cat", "Antibiotic"]).size().reset_index(name="COG_cat_freq")
data_table_overall = pd.pivot_table(data_grouped_overall, values='COG_cat_freq', index='Antibiotic', columns='Overall_cat', fill_value=0)
data_cols_overall = sorted(data_grouped_overall["Overall_cat"].unique())
plt.style.use('seaborn-pastel')
data_table_overall[data_cols_overall] = data_table_overall[data_cols_overall].div(data_table_overall[data_cols_overall].sum(axis=1), axis=0).multiply(100)
data_table_overall.plot(kind='barh', stacked=True, legend=True, width=1.0)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('Stacked_barchart_3_cat.pdf')

# subset into 3 main categories:
Cellular_processes = data_reduced[data_reduced["Overall_cat"] == "Cellular_processess_and_signaling"]
Information_storage = data_reduced[data_reduced["Overall_cat"] == "Information_storage_and_processing"]
Metabolism = data_reduced[data_reduced["Overall_cat"] == "Metabolism"]

# Cellular processes stacked barchart:
Cellular_processes_grouped = Cellular_processes.groupby(["COG_cat", "Antibiotic"]).size().reset_index(name="COG_cat_freq")
Cellular_processes_table = pd.pivot_table(Cellular_processes_grouped, values='COG_cat_freq', index='Antibiotic', columns='COG_cat', fill_value=0)
Cellular_processes_cols = sorted(Cellular_processes_grouped["COG_cat"].unique())
plt.style.use('seaborn-poster')
Cellular_processes_table[Cellular_processes_cols] = Cellular_processes_table[Cellular_processes_cols].div(Cellular_processes_table[Cellular_processes_cols].sum(axis=1), axis=0).multiply(100)
Cellular_processes_table.plot(kind='barh', stacked=True, edgecolor = 'white', legend=True, width=1.0, cmap='Pastel2', title='Cellular processes')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('Stacked_barchart_cellular_processess.pdf')

# Information storage stacked barchart:
Information_storage_grouped = Information_storage.groupby(["COG_cat", "Antibiotic"]).size().reset_index(name="COG_cat_freq")
Information_storage_table = pd.pivot_table(Information_storage_grouped, values='COG_cat_freq', index='Antibiotic', columns='COG_cat', fill_value=0)
Information_storage_cols = sorted(Information_storage_grouped["COG_cat"].unique())
plt.style.use('seaborn-poster')
Information_storage_table[Information_storage_cols] = Information_storage_table[Information_storage_cols].div(Information_storage_table[Information_storage_cols].sum(axis=1), axis=0).multiply(100)
Information_storage_table.plot(kind='barh', stacked=True, edgecolor = 'white', legend=True, width=1.0, cmap='Set3_r', title='Information storage')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('Stacked_barchart_information_storage.pdf')

# Metabolism stacked barchart:
Metabolism_grouped = Metabolism.groupby(["COG_cat", "Antibiotic"]).size().reset_index(name="COG_cat_freq")
Metabolism_table = pd.pivot_table(Metabolism_grouped, values='COG_cat_freq', index='Antibiotic', columns='COG_cat', fill_value=0)
Metabolism_cols = sorted(Metabolism_grouped["COG_cat"].unique())
plt.style.use('seaborn-poster')
Metabolism_table[Metabolism_cols] = Metabolism_table[Metabolism_cols].div(Metabolism_table[Metabolism_cols].sum(axis=1), axis=0).multiply(100)
Metabolism_table.plot(kind='barh', stacked=True, legend=True, edgecolor = 'white', width=1.0, cmap='Pastel1', title='Metabolism')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('Stacked_barchart_metabolism.pdf')