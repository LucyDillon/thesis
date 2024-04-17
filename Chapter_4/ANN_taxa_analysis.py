import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test data file ')
    parser.add_argument('-test_file', action='store', dest='test_file', required=True,
                            help='test data file')
    parser.add_argument('-training_data', action='store', dest='training_data',
                            help='training data', required=True)
    parser.add_argument('-result_file', action='store', dest='result_file',
                            help='result file', required=True)

    options = parser.parse_args()

    # DATA must be in a table format i.e. one genome per line with the phenotype:
    data = pd.read_csv(options.training_data)

    # drop first column
    data = data.drop(columns=["Unnamed: 0", "antibiotic", "Genus"])

    # Convert values to binary (0 or 1)
    binary_threshold = 1
    binary_columns = data.columns[data.columns.get_loc('genome_id')+1:data.columns.get_loc('add_phenotype')]
    data[binary_columns] = data[binary_columns].apply(lambda x: (x >= binary_threshold).astype(int))

    # Separate features (X) and target labels (y)
    X = data.drop(['genome_id', 'add_phenotype'], axis=1)  # Exclude non-feature columns
    y = LabelEncoder().fit_transform(data['add_phenotype'])

    # Split the data into 90% training and 10% testing
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

    # Create a new model
    model = Sequential()
    model.add(Dense(64, activation='relu', input_dim=X.shape[1]))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))  # 1 output neuron for binary classification

    # Compile the model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Train the model
    model.fit(X_train, y_train, epochs=10, batch_size=32)

    # Evaluate the model on the test set
    y_pred_prob = model.predict(X_test)
    y_pred = (y_pred_prob > 0.5).astype(int)

    # Calculate and print metrics
    test_accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    conf_matrix = confusion_matrix(y_test, y_pred)

    print(f"Test accuracy: {test_accuracy}")
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"Confusion matrix:\n{conf_matrix}")

    # Save the model
    model.save(f"{options.result_file}.ANN.h5")

    # Read in test data file
    test_data = pd.read_csv(options.test_file)

    # drop the first column
    test_data = test_data.drop(columns=["Unnamed: 0", "antibiotic", "Genus"])

    binary_threshold = 1
    binary_columns = test_data.columns[test_data.columns.get_loc('genome_id')+1:test_data.columns.get_loc('add_phenotype')]
    test_data[binary_columns] = test_data[binary_columns].apply(lambda x: (x >= binary_threshold).astype(int))

    # Extract genome_id
    genome_ids = test_data['genome_id'].values

    # Separate features (X) and target labels (y)
    X_test = test_data.drop(['genome_id', 'add_phenotype'], axis=1)
    y_test = LabelEncoder().fit_transform(test_data['add_phenotype'])

    # Ensure that X_test is reshaped appropriately
    X_test = X_test.values.reshape(X_test.shape[0], X_test.shape[1])

    # Load the trained model
    model = tf.keras.models.load_model(f"{options.result_file}.ANN.h5")

    # Make predictions on the test data
    predictions = model.predict(X_test)
    predicted_classes = (predictions > 0.5).astype(int)  # Adjust the threshold as needed
    phenotype_mapping = {0: 'Resistant', 1: 'Susceptible'} 
    # Create a DataFrame with genome_id and predictions
    result_df = pd.DataFrame({'genome_id': genome_ids, 'prediction': predicted_classes.flatten()})
    result_df['predicted_phenotype'] = result_df['prediction'].map(phenotype_mapping)

    # Save the result predictions:
    result_df.to_csv(f"{options.result_file}.ANN.csv", index=False)
    
