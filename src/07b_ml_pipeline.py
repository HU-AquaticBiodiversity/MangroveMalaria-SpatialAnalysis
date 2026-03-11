# ==============================================================================
# Script Name: Machine Learning Pipeline for Binary Classification
# Description: Trains and evaluates Logistic Regression, XGBoost, and a Keras 
#              Neural Network using GridSearchCV. Includes data preprocessing,
#              feature scaling, and performance evaluation via ROC AUC.
# ==============================================================================

## -----------------------------------------------------------------------------
## 1. LOAD LIBRARIES & CONFIGURE ENVIRONMENT
## -----------------------------------------------------------------------------
import os
# Disable TensorFlow warning logs and oneDNN optimizations for cleaner console output
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import tensorflow as tf
import pandas as pd
import numpy as np
import shap

# Scikit-Learn Preprocessing & Evaluation
from sklearn import preprocessing
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import roc_auc_score, classification_report

# Machine Learning Models
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from scikeras.wrappers import KerasClassifier
from keras.models import Sequential
from keras.layers import Input, Dense, Reshape

## -----------------------------------------------------------------------------
## 2. IMPORT DATA & PREPROCESSING
## -----------------------------------------------------------------------------
df = pd.read_csv('./data/ML_input_file.csv')



# ONE-HOT ENCODING:
# Machine learning models require numeric input. We convert the categorical 
# 'method' variable into distinct binary columns (0s and 1s).
encoder = OneHotEncoder(sparse_output=False).set_output(transform="pandas")
method_encoded = encoder.fit_transform(df[['method']])

# Append the newly encoded columns and drop the original text-based column
df = pd.concat([df, method_encoded], axis=1).drop(columns=['method'])

# GEOMETRIC TRANSFORMATION:
# Convert latitude and longitude to radians. This is highly recommended if your 
# downstream spatial algorithms calculate distances (e.g., Haversine formula).
df['lat'] = np.radians(df['lat'])
df['lon'] = np.radians(df['lon']) # FIX: Corrected duplicate 'lat' conversion

print("Data successfully loaded and preprocessed:")
print(df.head())

## -----------------------------------------------------------------------------
## 3. DEFINE FEATURES (X) AND TARGET (y)
## -----------------------------------------------------------------------------
# Isolate the predictor variables by dropping the target and survey weights
X = df.drop(columns=['infected', 'weight'])

# FEATURE SCALING:
# Standardize features (mean=0, variance=1) so that variables with large ranges
# do not disproportionately dominate the models (especially Logistic Reg & Keras).
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = pd.DataFrame(scaler.transform(X), columns=X.columns)

# Re-attach the sample weights temporarily for the train/test split
X_scaled['weight'] = df['weight'].to_list()

# Extract target variable as a 1D numpy array
y = df[['infected']].values.ravel()
    
# TRAIN/TEST SPLIT:
# Hold out 20% of the dataset to evaluate the final models on unseen data.
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, random_state=123
)

# Set up SHAP masker for potential feature importance extraction later
masker = shap.maskers.Independent(data=X_train.drop(columns=['weight']))

# Separate the survey sample weights from the training/testing features
weight = X_train['weight']

print("\nTraining features subset:")
print(X_train.head())

## -----------------------------------------------------------------------------
## 4. MODEL 1: LOGISTIC REGRESSION
## -----------------------------------------------------------------------------
print("\n--- Training Logistic Regression ---")
model_lr = LogisticRegression(random_state=0, max_iter=10000)

# HYPERPARAMETER GRID:
# 'C' controls regularization strength (smaller = stronger penalty on complexity).
param_lr = [{
    'C': np.logspace(-4, 4, 20),
    'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga']
}]

# Initialize and run Grid Search with 5-fold cross-validation
lr = GridSearchCV(model_lr, 
                  param_grid=param_lr, 
                  scoring='roc_auc',
                  cv=5, 
                  n_jobs=-1)

lr.fit(X_train.drop(columns=['weight']), y_train, sample_weight=weight)

print("Best Parameters:", lr.best_estimator_)
print("Best CV Score (ROC AUC):", lr.best_score_)

# Evaluate on the held-out test set

lr_pred = lr.predict(X_test.drop(columns=['weight']))
print("Test ROC AUC:", roc_auc_score(y_test, lr_pred))

## -----------------------------------------------------------------------------
## 5. MODEL 2: XGBOOST
## -----------------------------------------------------------------------------
print("\n--- Training XGBoost ---")
model_xgb = XGBClassifier(learning_rate=0.02, n_estimators=600, 
                          objective='binary:logistic', verbosity=0)

# HYPERPARAMETER GRID:
# Tuning tree constraints and data subsampling to prevent overfitting.
params_xgb = {
    'min_child_weight': [1, 5, 10],
    'gamma': [0.5, 1, 1.5, 2, 5],
    'subsample': [0.6, 0.8, 1.0],
    'colsample_bytree': [0.6, 0.8, 1.0],
    'max_depth': [3, 4, 5]
}

xgb = GridSearchCV(estimator=model_xgb, param_grid=params_xgb, 
                   scoring='roc_auc', n_jobs=-1, cv=5)

xgb.fit(X_train.drop(columns=['weight']), y_train, sample_weight=weight)

print("Best Parameters:", xgb.best_estimator_)
print("Best CV Score (ROC AUC):", xgb.best_score_)

# Evaluate on the held-out test set
xgb_pred = xgb.predict(X_test.drop(columns=['weight']))
print("Test ROC AUC:", roc_auc_score(y_test, xgb_pred))
print("Classification Report:\n", classification_report(xgb_pred, y_test))

## -----------------------------------------------------------------------------
## 6. MODEL 3: KERAS NEURAL NETWORK
## -----------------------------------------------------------------------------
print("\n--- Training Keras Neural Network ---")



[Image of Artificial Neural Network architecture]


# Define the architecture of the neural network dynamically
def model_build_fn(unit):
     model = Sequential([
         # Input layer expecting 21 numeric features
         Input(shape=(21,)),
         
         # Hidden layer
         Dense(units=unit, activation="relu", name="layer1"),
         
         # Output layer: 
         # FIX: Changed activation from "relu" to "sigmoid". 
         # Binary crossentropy requires outputs to be probabilities (0 to 1).
         Dense(1, activation="sigmoid", name="layer2"),
     ])
     return model
 
# Wrap the Keras model so it behaves like a Scikit-Learn estimator
keras_clf = KerasClassifier(model=model_build_fn, 
                            optimizer='adam',
                            loss='binary_crossentropy',
                            metrics=['accuracy'],
                            unit='none',  # Placeholder to be tuned by GridSearchCV
                            verbose=0)

# HYPERPARAMETER GRID:
# Tuning batch sizes, hidden layer node count (units), and training iterations (epochs)
param_grid = {
   'batch_size': [100, 20, 50, 25, 32], 
   'unit': [5, 6, 10, 11, 12, 15],
   'epochs': [100, 200, 300, 400],
}

grid = GridSearchCV(estimator=keras_clf,
                    param_grid=param_grid,
                    scoring='roc_auc', 
                    cv=5, 
                    n_jobs=-1,
                    verbose=1)

# Fit the grid search to the training data
keras_results = grid.fit(X_train.drop(columns=['weight']), 
                         y_train,
                         sample_weight=weight)

print("Best Parameters:", keras_results.best_params_)
print("Best CV Score (ROC AUC):", keras_results.best_score_)

# Evaluate on the held-out test set
keras_pred = keras_results.predict(X_test.drop(columns=['weight']))
print("Test ROC AUC:", roc_auc_score(y_test, keras_pred))
