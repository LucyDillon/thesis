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

    # Perform a single train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)  # You can adjust the test_size

    # Reshape data for 1D CNN (assuming each row is a sequence)
    X_train = X_train.values.reshape(X_train.shape[0], X_train.shape[1], 1)
    X_test = X_test.values.reshape(X_test.shape[0], X_test.shape[1], 1)

    # Create a new 1D CNN model
    model = Sequential()
    model.add(Conv1D(32, kernel_size=3, activation='relu', input_shape=(X_train.shape[1], 1)))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))

    # Compile the 1D CNN model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Train the 1D CNN model
    model.fit(X_train, y_train, epochs=10, batch_size=32)

    # Evaluate the 1D CNN model on the test set
    test_loss, test_accuracy = model.evaluate(X_test, y_test)

    # Calculate precision and recall
    y_pred_prob = model.predict(X_test)
    y_pred = (y_pred_prob > 0.5).astype(int)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)

    # Calculate and print the confusion matrix
    conf_matrix = confusion_matrix(y_test, y_pred)

    # Print results
    print(f"Test accuracy: {test_accuracy}")
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"Confusion matrix:\n{conf_matrix}")

    # Save the best model outside the loop
    if model is not None:
        model.save(f"{options.result_file}.h5")
        print(f"Best model saved with test accuracy: {test_accuracy}")
    else:
        print("No model saved as no fold was processed.")

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
    model = tf.keras.models.load_model(f"{options.result_file}.h5")

    # Make predictions on the test data
    predictions = model.predict(X_test)
    predicted_classes = (predictions > 0.5).astype(int)  # Adjust the threshold as needed
    phenotype_mapping = {0: 'Resistant', 1: 'Susceptible'} 
    # Create a DataFrame with genome_id and predictions
    result_df = pd.DataFrame({'genome_id': genome_ids, 'prediction': predicted_classes.flatten()})
    result_df['predicted_phenotype'] = result_df['prediction'].map(phenotype_mapping)

    # Save the result predictions:
    result_df.to_csv(options.result_file, index=False)
    