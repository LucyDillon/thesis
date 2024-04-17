import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score


# DATA must be in a table format i.e. one genome per line with the phenotype:
data = pd.read_csv('ampicillin.csv')
# Convert values to binary (0 or 1)
binary_threshold = 1
binary_columns = data.columns[data.columns.get_loc('genome_id')+1:data.columns.get_loc('add_phenotype')]
data[binary_columns] = data[binary_columns].apply(lambda x: (x >= binary_threshold).astype(int))

# Separate features (X) and target labels (y)
X = data.drop(['genome_id', 'add_phenotype'], axis=1)  # Exclude non-feature columns
y = LabelEncoder().fit_transform(data['add_phenotype'])

# Number of folds for cross-validation
num_folds = 10
kf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)

accuracies = []
precisions = []
recalls = []
conf_matrices = []

for train_idx, test_idx in kf.split(X, y):
    X_train_fold, X_test_fold = X.iloc[train_idx], X.iloc[test_idx]
    y_train_fold, y_test_fold = y[train_idx], y[test_idx]

    # Create a new model for each fold (optional)
    model = Sequential()
    model.add(Dense(64, activation='relu', input_dim=X.shape[1]))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))  # 1 output neuron for binary classification

    # Compile the model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Train the model on the current fold
    model.fit(X_train_fold, y_train_fold, epochs=10, batch_size=32)

    # Evaluate the model on the test fold
    y_pred_prob = model.predict(X_test_fold)
    y_pred = (y_pred_prob > 0.5).astype(int)

    test_loss, test_accuracy = model.evaluate(X_test_fold, y_test_fold)
    accuracies.append(test_accuracy)

    # Calculate precision and recall
    precision = precision_score(y_test_fold, y_pred)
    recall = recall_score(y_test_fold, y_pred)
    precisions.append(precision)
    recalls.append(recall)

    # Calculate and print the confusion matrix
    conf_matrix = confusion_matrix(y_test_fold, y_pred)
    conf_matrices.append(conf_matrix)

    print(f"Test accuracy for fold: {test_accuracy}")
    print(f"Precision for fold: {precision}")
    print(f"Recall for fold: {recall}")
    print(f"Confusion matrix for fold:\n{conf_matrix}")

# Calculate and print the average accuracy, precision, and recall across all folds
average_accuracy = np.mean(accuracies)
average_precision = np.mean(precisions)
average_recall = np.mean(recalls)

print(f"Average Test Accuracy: {average_accuracy}")
print(f"Average Precision: {average_precision}")
print(f"Average Recall: {average_recall}")

# Save the model
model.save('ampicillin_egg_model_K_fold_bin.h5')
