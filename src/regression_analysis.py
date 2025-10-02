## load libraries ##
import tensorflow as tf
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier
from scikeras.wrappers import KerasClassifier
from sklearn.metrics import classification_report
from keras.models import Sequential
from keras.layers import Dense
import pandas as pd
import numpy as np
import shap

##---------------------------------------------------##
## import data frame                                 ##
##---------------------------------------------------##
df = pd.read_csv('./data/ML_input_file.csv')

# transform categorical variable
encoder = OneHotEncoder(sparse_output=False).set_output(transform="pandas")
method_encoded = encoder.fit_transform(df[['method']])
df = pd.concat([df, method_encoded],
               axis = 1).drop(columns = ['method'])

# transform degrees to radian for better performance
df['lat'] = np.radians(df['lat'])
df['lat'] = np.radians(df['lat'])
print(df.head())

##---------------------------------------------------##
## defining feature matrix(X) and response vector(y) ##
##---------------------------------------------------##
    
X = df.drop(columns = ['infected', 'weight'])
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = pd.DataFrame(scaler.transform(X), columns = X.columns)
X_scaled['weight'] = df['weight'].to_list()

# needs to be transformed to 1d array
y = df[['infected']].values.ravel()
    
# split between training and test data
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y,test_size=0.2, random_state=123
        )
masker = shap.maskers.Independent(data =
                                  X_train.drop(columns = ['weight']))
weight = X_train['weight']

##---------------------------------------------------##
## Logistic Regression                               ##
##---------------------------------------------------##
    
# Starting model
model_lr = LogisticRegression(random_state=0, max_iter = 1000)
# make parameter grid
param_lr = [
        {
            'C' : np.logspace(-4,4,20),
            'solver': ['lbfgs','newton-cg','liblinear','sag','saga']
            }
        ]

# define grid search
lr = GridSearchCV(model_lr, 
                  param_grid = param_lr, 
                  scoring='roc_auc',
                  cv = 5, 
                  n_jobs = 10)
# run grid search
lr.fit(X_train.drop(columns = ['weight']), y_train, sample_weight = weight)
# show best parameter
print("The best parameters are:")
print(lr.best_estimator_)
print(lr.best_score_)

#optimised model
lr_pred = lr.predict(X_test.drop(columns = ['weight']))
       
# print classification report 
print("roc auc for test data:")
print(roc_auc_score(y_test, lr_pred)) 

# shap values
explainer_lr = shap.LinearExplainer(lr.best_estimator_, masker = masker)
shap_lr = explainer_lr(X_train.drop(columns = ['weight']))

# visualize the first prediction's explanation
shap.plots.waterfall(shap_lr[0], max_display = 44)
shap.summary_plot(shap_lr, X_train.drop(columns = ['weight']), 
                  plot_type= 'bar', max_display= 43)

##---------------------------------------------------##
## XG Boost ##
##---------------------------------------------------##
model_xgb = XGBClassifier(learning_rate=0.02, n_estimators=600, 
                    objective='binary:logistic',
                    verbosity=0)

# A parameter grid for XGBoost
params_xgb = {
        'min_child_weight': [1, 5, 10],
        'gamma': [0.5, 1, 1.5, 2, 5],
        'subsample': [0.6, 0.8, 1.0],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'max_depth': [3, 4, 5]
        }

xgb = GridSearchCV(estimator=model_xgb, param_grid=params_xgb, scoring='roc_auc',
                    n_jobs=10, cv= 5)

# run grid search
xgb.fit(X_train.drop(columns = ['weight']), y_train, sample_weight = weight)

# show best parameter
print("The best parameters are:")
print(xgb.best_estimator_)
print(xgb.best_score_)

#optimised model
xgb_pred = xgb.predict(X_test.drop(columns = ['weight']))
   
# print classification report 
print("roc auc for test data:")
print(roc_auc_score(y_test, xgb_pred))
print(classification_report(xgb_pred, y_test))

# shap values
#explainer_xgb = shap.TreeExplainer(xgb.best_estimator_)
#shap_xgb = explainer_xgb(X_train.drop(columns = ['weight']))

#shap.summary_plot(shap_xgb, X_train.drop(columns = ['weight']), plot_type='bar')

##---------------------------------------------------##
## Keras ##
##---------------------------------------------------##
# Define Keras model function

# Define a function to create the Keras model
def create_model(units=10, activation='relu'):
   model = Sequential()
   model.add(Dense(units=units, 
                   activation=activation,
                   input_shape = (33,)))
   model.add(Dense(units=1, activation='sigmoid'))
   model.compile(optimizer='adam', loss='binary_crossentropy', 
                 metrics=['accuracy'])
   return model

# Create the KerasClassifier object
model = KerasClassifier(build_fn=create_model)

# Define the hyperparameter grid to search over
param_grid = {
   'model__units': [5,6, 10, 11, 12, 15],
   'model__activation': ['relu', 'sigmoid'],
   'epochs':[200, 100, 300, 400],
}

# Create the GridSearchCV object
grid = GridSearchCV(estimator=model, param_grid = param_grid, 
                    scoring='roc_auc', cv=5, n_jobs = -1)

# Fit the GridSearchCV object to the training data
keras_results = grid.fit(X_train.drop(columns = ['weight']), y_train, 
             sample_weight = weight)

# Print the best parameters and score
print("Best Parameters: ", keras_results.best_params_)
print("Best Score: ", keras_results.best_score_)

#optimised model
keras_pred = keras_results.predict(X_test.drop(columns = ['weight']))
   
# print classification report 
print("roc auc for test data:")
print(roc_auc_score(y_test, xgb_pred))

# shap values
#explainer_xgb = shap.TreeExplainer(xgb.best_estimator_)
#shap_xgb = explainer_xgb(X_test.drop(columns = ['weight']))

# visualize the first prediction's explanation
#shap.plots.waterfall(shap_xgb[0], max_display = 30)

### additional code

for r in range(1,51):    

       
    df_select = df[df["buffer"] == r]
    
    X = df_select.drop(columns = ['infected', 'date', 'weight'])
    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = pd.DataFrame(scaler.transform(X), columns = X.columns)
    X_scaled['weight'] = df_select['weight'].to_list()

    # needs to be transformed to 1d array
    y = df_select[['infected']].values.ravel()
    
    # split between training and test data
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y,test_size=0.1, random_state=123
        )
    
    masker = shap.maskers.Independent(data = 
                                      X_train.drop(columns = ['weight']))
    weight = X_train['weight']
            
    ##---------------------------------------------------##
    ## Logistic Regression                               ##
    ##---------------------------------------------------##
    
    # Starting model
    model_lr = LogisticRegression(random_state=0, max_iter = 1000)
    # make parameter grid
    param_lr = [
        {
        'C' : np.logspace(-4,4,20)
        }
    ]

    # define grid search
    lr = GridSearchCV(model_lr, 
                       param_grid = param_lr, 
                       scoring='roc_auc',
                       cv = 5, 
                       n_jobs = 10)
    # run grid search
    lr.fit(X_train.drop(columns = ['weight']), y_train, sample_weight = weight)
    # show best parameter
    print("The best parameters are:")
    print(lr.best_estimator_)
    print(lr.best_score_)

    #optimised model
    lr_pred = lr.predict(X_test.drop(columns = ['weight']))
       
    # print classification report 
    print("roc auc for test data:")
    print(roc_auc_score(y_test, lr_pred)) 
