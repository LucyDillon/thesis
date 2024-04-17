import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# List of antibiotics
antibiotics = ['amikacin', 'amoxicillin', 'ampicillin', 'aztreonam', 'ceftriaxone', 'cefepime', 'chloramphenicol', 'ciprofloxacin', 'clindamycin', 'colistin', 'doripenem', 'ertapenem', 'erythromycin', 'fosfomycin', 'gentamicin', 'imipenem', 'levofloxacin', 'meropenem', 'moxifloxacin', 'nitrofurantoin', 'tetracycline', 'tigecycline', 'tobramycin']

# List of taxonomic levels
taxonomic_levels = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum']

# Loop over antibiotics
for antibiotic in antibiotics:
    # Read data from file
    file_name = f"{antibiotic}_tree_traversal.txt"
    data = pd.read_csv(file_name, sep='\t')

    # Loop over taxonomic levels
    for level in taxonomic_levels:
        # Group and calculate percentages
        data['Group'] = data.groupby(['Node_no', 'Node_label']).ngroup() + 1
        data_grouped = data.groupby('Group')[level].value_counts(normalize=True).unstack() * 100

        # Plotting
        plt.style.use('seaborn-colorblind')
        ax = data_grouped.plot(kind='bar', stacked=True, title=f"{antibiotic.capitalize()} pathways to resistance - Eggnog model by {level}", edgecolor='white', width=1.0)
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize="15")
        sns.set_style("ticks")
        plt.xticks(fontsize=20)
        ax.set_xlabel('Pathway', fontsize=20)
        ax.set_ylabel('Percentage (%)', fontsize=20)
        ax.set_title(f"{antibiotic.capitalize()} pathways to resistance - Eggnog model by {level}", fontsize=25)
        fig = ax.get_figure()
        fig.set_size_inches(22, 6)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        fig.savefig(f"Eggnog_{antibiotic.capitalize()}_{level}_pathways.svg", format='svg')
        plt.close(fig)  # Close the figure to avoid overlapping plots
