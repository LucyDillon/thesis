import os
import pandas as pd

# Specify the directory containing the pairwise files
directory = "./"

# Create an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

# Loop through each pairwise file in the directory
for filename in os.listdir(directory):
    if filename.endswith("_contrib_to_phenotype.txt") or filename.endswith("_presense_contrib_to_phenotype.txt"):
        file_path = os.path.join(directory, filename)
        
        # Read the current file
        current_df = pd.read_csv(file_path, sep='\t')
        
        # Extract the drug names from the filename
        drugs = filename.split('_')[2:4]
        drug_pair = '_'.join(drugs)
        
        # Extract presence or absence from the filename
        contribution_type = filename.split('_')[6]
        
        # Add new columns for drug pair and presence/absence in the dataframe
        current_df['Drug_pair'] = drug_pair
        current_df['Contribution_type'] = contribution_type
        
        # Append the current data to the merged DataFrame
        merged_df = pd.concat([merged_df, current_df], ignore_index=True)

# Remove rows where any of the values (SS, SR, RS, RR) is NaN
merged_df = merged_df.dropna(subset=['SS', 'SR', 'RS', 'RR'])

# Sort the merged DataFrame by Gene_name, Drug_pair, and Contribution_type
merged_df.sort_values(['Gene_name', 'Drug_pair', 'Contribution_type'], inplace=True)

# Print the result or save it to a new file
print(merged_df)
merged_df.to_csv("merged_file.csv", index=False)
