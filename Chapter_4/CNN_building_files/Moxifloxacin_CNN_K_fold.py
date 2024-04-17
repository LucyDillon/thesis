import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score

# DATA must be in a table format i.e. one genome per line with the phenotype:
data = pd.read_csv('moxifloxacin.csv')
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

all_accuracies = []
all_precisions = []
all_recalls = []
all_conf_matrices = []

best_accuracy = 0
best_model = None

for fold, (train_idx, test_idx) in enumerate(kf.split(X, y), 1):
    X_train_fold, X_test_fold = X.iloc[train_idx].values, X.iloc[test_idx].values
    y_train_fold, y_test_fold = y[train_idx], y[test_idx]

    # Reshape data for 1D CNN (assuming each row is a sequence)
    X_train_fold = X_train_fold.reshape(X_train_fold.shape[0], X_train_fold.shape[1], 1)
    X_test_fold = X_test_fold.reshape(X_test_fold.shape[0], X_test_fold.shape[1], 1)

    # Create a new 1D CNN model for each fold
    model = Sequential()
    model.add(Conv1D(32, kernel_size=3, activation='relu', input_shape=(X_train_fold.shape[1], 1)))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))

    # Compile the 1D CNN model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Train the 1D CNN model on the current fold
    model.fit(X_train_fold, y_train_fold, epochs=10, batch_size=32)

    # Evaluate the 1D CNN model on the test fold
    test_loss, test_accuracy = model.evaluate(X_test_fold, y_test_fold)

    # Save the model if it has the highest test accuracy
    if test_accuracy > best_accuracy:
        best_accuracy = test_accuracy
        best_model = model

    # Calculate precision and recall
    y_pred_prob = model.predict(X_test_fold)
    y_pred = (y_pred_prob > 0.5).astype(int)
    precision = precision_score(y_test_fold, y_pred)
    recall = recall_score(y_test_fold, y_pred)

    # Calculate and print the confusion matrix
    conf_matrix = confusion_matrix(y_test_fold, y_pred)

    # Print results for the current fold
    print(f"Fold {fold}: Test accuracy: {test_accuracy}")
    print(f"Fold {fold}: Precision: {precision}")
    print(f"Fold {fold}: Recall: {recall}")
    print(f"Fold {fold}: Confusion matrix:\n{conf_matrix}\n")

    # Save metrics for later averaging
    all_accuracies.append(test_accuracy)
    all_precisions.append(precision)
    all_recalls.append(recall)
    all_conf_matrices.append(conf_matrix)

# Calculate and print the average accuracy, precision, recall, and confusion matrix across all folds
average_accuracy = np.mean(all_accuracies)
average_precision = np.mean(all_precisions)
average_recall = np.mean(all_recalls)
average_conf_matrix = np.mean(all_conf_matrices, axis=0)

print(f"Average Test Accuracy: {average_accuracy}")
print(f"Average Precision: {average_precision}")
print(f"Average Recall: {average_recall}")
print(f"Average Confusion Matrix:\n{average_conf_matrix}")

# Save the best model outside the loop
if best_model is not None:
    best_model.save('moxifloxacin_CNN_K_fold.h5')
    print(f"Best model saved with test accuracy: {best_accuracy}")
else:
    print("No model saved as no fold was processed.")
