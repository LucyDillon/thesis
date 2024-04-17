import pandas as pd
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predictions data')
    parser.add_argument('-predictions_file', action='store', dest='predictions_file', required=True,
                            help='predictions file')
    parser.add_argument('-actual_pheno', action='store', dest='actual_pheno',
                            help='actual phenotype', required=True)
    parser.add_argument('-model_name', action='store', dest='model_name',
                            help='model name', required=True)
    
    options = parser.parse_args()

    predictions = pd.read_csv(options.predictions_file)
    Actual_pheno =pd.read_csv(options.actual_pheno, sep=' ')

    merged_data = pd.merge(left=predictions, right=Actual_pheno, on='genome_id')
    merged_data['agree'] = np.where(merged_data['predicted_phenotype'] == merged_data['add_phenotype'], 1, 0)

    percentage_agree = merged_data['agree'].mean() * 100
    print(f"{options.model_name}: {percentage_agree}% correct predi