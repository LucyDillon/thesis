#import packages:
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.model_selection import StratifiedKFold,train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, multilabel_confusion_matrix


AB_files = ["E_coli_amikacin_ampicillin_MDR.csv", "E_coli_ampicillin_cefepime_MDR.csv", "E_coli_aztreonam_doripenem_MDR.csv", "E_coli_cefepime_levofloxacin_MDR.csv", "E_coli_ciprofloxacin_doripenem_MDR.csv", "E_coli_doripenem_tigecycline_MDR.csv", "E_coli_imipenem_levofloxacin_MDR.csv",
"E_coli_amikacin_aztreonam_MDR.csv", "E_coli_ampicillin_ceftriaxone_MDR.csv", "E_coli_aztreonam_ertapenem_MDR.csv", "E_coli_cefepime_meropenem_MDR.csv", "E_coli_ciprofloxacin_ertapenem_MDR.csv", "E_coli_doripenem_tobramycin_MDR.csv", "E_coli_imipenem_meropenem_MDR.csv",
"E_coli_amikacin_cefepime_MDR.csv", "E_coli_ampicillin_ciprofloxacin_MDR.csv", "E_coli_aztreonam_gentamicin_MDR.csv", "E_coli_cefepime_nitrofurantoin_MDR.csv", "E_coli_ciprofloxacin_gentamicin_MDR.csv", "E_coli_ertapenem_gentamicin_MDR.csv", "E_coli_imipenem_tigecycline_MDR.csv",
"E_coli_amikacin_ceftriaxone_MDR.csv", "E_coli_ampicillin_doripenem_MDR.csv", "E_coli_aztreonam_imipenem_MDR.csv", "E_coli_cefepime_tigecycline_MDR.csv", "E_coli_ciprofloxacin_imipenem_MDR.csv", "E_coli_ertapenem_imipenem_MDR.csv", "E_coli_imipenem_tobramycin_MDR.csv",
"E_coli_amikacin_ciprofloxacin_MDR.csv", "E_coli_ampicillin_ertapenem_MDR.csv", "E_coli_aztreonam_levofloxacin_MDR.csv", "E_coli_cefepime_tobramycin_MDR.csv", "E_coli_ciprofloxacin_levofloxacin_MDR.csv",	"E_coli_ertapenem_levofloxacin_MDR.csv", "E_coli_levofloxacin_meropenem_MDR.csv",
"E_coli_amikacin_doripenem_MDR.csv", "E_coli_ampicillin_gentamicin_MDR.csv", "E_coli_aztreonam_meropenem_MDR.csv", "E_coli_ceftriaxone_ciprofloxacin_MDR.csv",	"E_coli_ciprofloxacin_meropenem_MDR.csv", "E_coli_ertapenem_meropenem_MDR.csv", "E_coli_levofloxacin_nitrofurantoin_MDR.csv",
"E_coli_amikacin_ertapenem_MDR.csv", "E_coli_ampicillin_imipenem_MDR.csv", "E_coli_aztreonam_nitrofurantoin_MDR.csv", "E_coli_ceftriaxone_doripenem_MDR.csv", "E_coli_ciprofloxacin_nitrofurantoin_MDR.csv",	"E_coli_ertapenem_nitrofurantoin_MDR.csv", "E_coli_levofloxacin_tigecycline_MDR.csv",
"E_coli_amikacin_gentamicin_MDR.csv", "E_coli_ampicillin_levofloxacin_MDR.csv", "E_coli_aztreonam_tigecycline_MDR.csv", "E_coli_ceftriaxone_ertapenem_MDR.csv", "E_coli_ciprofloxacin_tigecycline_MDR.csv",	"E_coli_ertapenem_tigecycline_MDR.csv", "E_coli_levofloxacin_tobramycin_MDR.csv",
"E_coli_amikacin_imipenem_MDR.csv", "E_coli_ampicillin_meropenem_MDR.csv", "E_coli_aztreonam_tobramycin_MDR.csv", "E_coli_ceftriaxone_gentamicin_MDR.csv", "E_coli_ciprofloxacin_tobramycin_MDR.csv", "E_coli_ertapenem_tobramycin_MDR.csv", "E_coli_meropenem_nitrofurantoin_MDR.csv",
"E_coli_amikacin_levofloxacin_MDR.csv", "E_coli_ampicillin_nitrofurantoin_MDR.csv",	"E_coli_cefepime_ceftriaxone_MDR.csv", "E_coli_ceftriaxone_imipenem_MDR.csv", "E_coli_doripenem_ertapenem_MDR.csv", "E_coli_gentamicin_imipenem_MDR.csv", "E_coli_meropenem_tigecycline_MDR.csv",
"E_coli_amikacin_meropenem_MDR.csv", "E_coli_ampicillin_tigecycline_MDR.csv", "E_coli_cefepime_ciprofloxacin_MDR.csv", "E_coli_ceftriaxone_levofloxacin_MDR.csv", "E_coli_doripenem_gentamicin_MDR.csv", "E_coli_gentamicin_levofloxacin_MDR.csv", "E_coli_meropenem_tobramycin_MDR.csv",
"E_coli_amikacin_nitrofurantoin_MDR.csv", "E_coli_ampicillin_tobramycin_MDR.csv", "E_coli_cefepime_doripenem_MDR.csv", "E_coli_ceftriaxone_meropenem_MDR.csv", "E_coli_doripenem_imipenem_MDR.csv", "E_coli_gentamicin_meropenem_MDR.csv", "E_coli_nitrofurantoin_tigecycline_MDR.csv",
"E_coli_amikacin_tigecycline_MDR.csv", "E_coli_aztreonam_cefepime_MDR.csv", "E_coli_cefepime_ertapenem_MDR.csv", "E_coli_ceftriaxone_nitrofurantoin_MDR.csv",	"E_coli_doripenem_levofloxacin_MDR.csv", "E_coli_gentamicin_nitrofurantoin_MDR.csv",	"E_coli_nitrofurantoin_tobramycin_MDR.csv",
"E_coli_amikacin_tobramycin_MDR.csv", "E_coli_aztreonam_ceftriaxone_MDR.csv", "E_coli_cefepime_gentamicin_MDR.csv", "E_coli_ceftriaxone_tigecycline_MDR.csv", "E_coli_doripenem_meropenem_MDR.csv", "E_coli_gentamicin_tigecycline_MDR.csv", "E_coli_tigecycline_tobramycin_MDR.csv",
"E_coli_ampicillin_aztreonam_MDR.csv", "E_coli_aztreonam_ciprofloxacin_MDR.csv", "E_coli_cefepime_imipenem_MDR.csv", "E_coli_ceftriaxone_tobramycin_MDR.csv", "E_coli_doripenem_nitrofurantoin_MDR.csv", "E_coli_gentamicin_tobramycin_MDR.csv"]


for j in AB_files:
    data = pd.read_csv(j)

    # Convert values to one-hot encoding for more than two phenotypes
    label_encoder = LabelEncoder()
    data['phenotype_encoded'] = label_encoder.fit_transform(data['phenotype'])

    one_hot_encoder = OneHotEncoder(sparse=False)
    phenotype_encoded = data['phenotype_encoded'].values.reshape(-1, 1)
    one_hot_encoded = one_hot_encoder.fit_transform(phenotype_encoded)
    one_hot_columns = [f'phenotype_{label}' for label in label_encoder.classes_]
    data[one_hot_columns] = pd.DataFrame(one_hot_encoded, columns=one_hot_columns)

    # Separate features (X) and target labels (y)
    X = data.drop(['genome_id', 'phenotype', 'phenotype_encoded'] + one_hot_columns, axis=1)
    y = data[one_hot_columns]

    # Perform a single train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)  # You can adjust the test_size

    # Reshape data for 1D CNN (assuming each row is a sequence)
    X_train = X_train.values.reshape(X_train.shape[0], X_train.shape[1], 1)
    X_test = X_test.values.reshape(X_test.shape[0], X_test.shape[1], 1)

        # Update the number of units in the output layer of your neural network
        # to match the number of unique phenotypes in your dataset
    output_units = len(label_encoder.classes_)

        # Create a new 1D CNN model
    model = Sequential()
    model.add(Conv1D(32, kernel_size=3, activation='relu', input_shape=(X_train.shape[1], 1)))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dense(output_units, activation='softmax'))  # Use softmax for multi-class classification
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

        # Train the 1D CNN model
    model.fit(X_train, y_train, epochs=10, batch_size=32)
        # Evaluate the 1D CNN model on the test set
    test_loss, test_accuracy = model.evaluate(X_test, y_test)

    # Save the best model outside the loop
    if model is not None:
        model.save(f"{j}_test_CNN.h5")
        print(f"Best model for {j} saved with test accuracy: {test_accuracy}")
    else:
        print(f"No model saved for {j} as no fold was processed.")

    y_pred = model.predict(X_test)

    # Convert the predicted probabilities to class labels
    y_pred_labels = np.argmax(y_pred, axis=1)
    y_test_labels = np.argmax(y_test.values, axis=1)

    # Calculate the confusion matrix
    conf_matrix = confusion_matrix(y_test_labels, y_pred_labels)

    # Print the confusion matrix as a table with labels
    print("Confusion matrix:")
    print("\t\t\tPredicted")
    print("\t\t\t" + "\t".join(label_encoder.classes_))
    for i, row in enumerate(conf_matrix):
        print(f"Actual {label_encoder.classes_[i]}\t", "\t".join(map(str, row)))




