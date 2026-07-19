data0='FullW4_piHRSD_selected'
traindata='NHRI'
print("traindata:", traindata)

import pandas as pd
from sklearn.model_selection import train_test_split
ds = pd.read_csv('/home/u7080475/SSRI/'+data0+'.csv')

test=1
X_train=ds.loc[ds['Site']==traindata]
X_test=ds.loc[ds['Site']!=traindata]
y_train=X_train['piHRSDW4']
y_test=X_test['piHRSDW4'] 

X_train=X_train.drop(['ID','Site','piHRSDW4','as.factor(medication)fluoxetine','as.factor(medication)citalopram','as.factor(medication)paroxetine','as.factor(medication)escitalopram'], axis=1, errors='ignore')
X_test=X_test.drop(['ID','Site','piHRSDW4','as.factor(medication)fluoxetine','as.factor(medication)citalopram','as.factor(medication)paroxetine','as.factor(medication)escitalopram'], axis=1, errors='ignore')


from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn import model_selection
from sklearn.model_selection import GridSearchCV, RepeatedKFold # 改用 RepeatedKFold
cv = model_selection.RepeatedKFold(n_splits=5, n_repeats=5, random_state=0) 

# ElasticNet 
from sklearn.linear_model import ElasticNet # ElasticNet
reg_model = ElasticNet()
param_grid = {'l1_ratio': [i * 0.1 for i in range(6,11)],'alpha': [10 ** i for i in range(-5, 3)], 'random_state': [0]}

if test==1: 
    reg_cv = GridSearchCV(reg_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 
else:
    reg_cv = GridSearchCV(reg_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 

reg_cv.fit(X_train, y_train)

print("Best Hyperparameters:", reg_cv.best_params_)
print("Best Score:", -reg_cv.best_score_)


#KNN
from sklearn.neighbors import KNeighborsRegressor
knn = KNeighborsRegressor()

param_grid = {'n_neighbors': [3, 5, 10, 15, 20, 25, 30, 40, 50, 60], 'weights': ['uniform', 'distance'], 'p': [1, 2]}
if test==1:
    knn_cv = GridSearchCV(knn, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 
else:
    knn_cv = GridSearchCV(knn, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 

knn_cv.fit(X_train, y_train)

print("Best Hyperparameters:", knn_cv.best_params_)
print("Best Score (MSE):", -knn_cv.best_score_)

# SVR
from sklearn.svm import SVR
svm = SVR()
param_grid = {'C': [10 ** i for i in range(2, -6, -1)], 'kernel': ['linear'], 'kernel': ['linear', 'rbf'], 'epsilon': [0.01, 0.1, 0.2, 0.5],'gamma': ['scale', 'auto']}
# param_grid = {'C': [10 ** i for i in range(0, -6, -1)], 'kernel': ['linear']}

if test==1:
    svm_cv = GridSearchCV(svm, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 
else:
    svm_cv = GridSearchCV(svm, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 

svm_cv.fit(X_train, y_train)

print("Best Hyperparameters:", svm_cv.best_params_)
print("Best Score (MSE):", -svm_cv.best_score_)
svmScore = -svm_cv.best_score_

# xgboost/XGBRegressor
from xgboost import XGBRegressor

xgboost_model = XGBRegressor(nthread=1) 
if test==1:
    param_grid = {
            'n_estimators': [100], 
            'booster': ['gbtree'], 
            'learning_rate': [0.1], 
            'gamma': [0], 
            'reg_alpha': [0], 
            'reg_lambda': [0], 
            'objective': ['reg:squarederror'],
            'colsample_bylevel': [0.7,0.9], 
            'colsample_bynode': [0.1,0.3], 
            'colsample_bytree': [0.1,0.3], 
            'subsample': [0.8],
            'min_child_weight': [5],  
            'max_depth': [1], 
            'random_state': [0]
            }
    xgboost_cv = GridSearchCV(xgboost_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 
else:
    param_grid = {
            'n_estimators': [500], 
            'booster': ['gbtree'], 
            'learning_rate': [0.1], 
            'gamma': [3], 
            'reg_alpha': [0], 
            'reg_lambda': [0],
            'objective': ['reg:squarederror'], 
            'colsample_bylevel': [0.1,0.3,0.5,0.7,0.9], 
            'colsample_bynode': [0.1,0.3,0.5,0.7,0.9], 
            'colsample_bytree': [0.1,0.3,0.5,0.7,0.9],
            'subsample': [0.6],
            'min_child_weight': [5],  
            'max_depth': [1,2,3,4,5], 
            'random_state': [0]
            }
    xgboost_cv = GridSearchCV(xgboost_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40) 

xgboost_cv.fit(X_train, y_train) 

print("Best Hyperparameters:", xgboost_cv.best_params_)
print("Best Score (MSE):", -xgboost_cv.best_score_)


#ensemble
from sklearn.ensemble import VotingRegressor
ensem = [('elasticnet', reg_cv.best_estimator_),('knn', knn_cv.best_estimator_),('svm', svm_cv.best_estimator_), ('xgboost', xgboost_cv.best_estimator_)]
ensemble = VotingRegressor(ensem)
ensemble.fit(X_train, y_train)
print("Test R2 Score:", ensemble.score(X_test, y_test))

import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
res = pd.DataFrame([
    ['ElasticNet', -reg_cv.best_score_,
     mean_squared_error(y_train, reg_cv.predict(X_train)),
     r2_score(y_train, reg_cv.predict(X_train)),
     mean_squared_error(y_test, reg_cv.predict(X_test)),
     r2_score(y_test, reg_cv.predict(X_test))],
    
    ['KNN', -knn_cv.best_score_,
     mean_squared_error(y_train, knn_cv.predict(X_train)),
     r2_score(y_train, knn_cv.predict(X_train)),
     mean_squared_error(y_test, knn_cv.predict(X_test)),
     r2_score(y_test, knn_cv.predict(X_test))],
    
    ['SVR', -svm_cv.best_score_,
     mean_squared_error(y_train, svm_cv.predict(X_train)),
     r2_score(y_train, svm_cv.predict(X_train)),
     mean_squared_error(y_test, svm_cv.predict(X_test)),
     r2_score(y_test, svm_cv.predict(X_test))],
    
    ['XGBoost', -xgboost_cv.best_score_,
     mean_squared_error(y_train, xgboost_cv.predict(X_train)),
     r2_score(y_train, xgboost_cv.predict(X_train)),
     mean_squared_error(y_test, xgboost_cv.predict(X_test)),
     r2_score(y_test, xgboost_cv.predict(X_test))],
    
    ['Ensemble', 'NA',
     mean_squared_error(y_train, ensemble.predict(X_train)),
     r2_score(y_train, ensemble.predict(X_train)),
     mean_squared_error(y_test, ensemble.predict(X_test)),
     r2_score(y_test, ensemble.predict(X_test))]
], columns=['Method', 'CV_MSE', 'MSE_train', 'R2_train', 'MSE_test', 'R2_test'])

print(res)
res.to_csv('/home/u7080475/SSRI/traditional_model/results/res'+data0+traindata+'.csv', index=False)


import pickle
save_dir = '/home/u7080475/SSRI/traditional_model/models/'
pickle.dump(reg_cv, open(save_dir + data0 + traindata + '_elasticnet.pkl', 'wb'))
pickle.dump(knn_cv, open(save_dir + data0 + traindata + '_knn.pkl', 'wb'))
pickle.dump(svm_cv, open(save_dir + data0 + traindata + '_svr.pkl', 'wb'))
pickle.dump(xgboost_cv, open(save_dir + data0 + traindata + '_xgboost.pkl', 'wb'))
pickle.dump(ensemble, open(save_dir + data0 + traindata + '_ensemble.pkl', 'wb'))


#========for SHAP
from multiprocessing import Pool
import pickle
import shap
import joblib

indx = res['R2_test'].idxmax()
best_method = res['Method'][indx]
print("best method(SHAP):", best_method)

file_suffix = best_method.lower()
loaded_model = joblib.load('/home/u7080475/SSRI/traditional_model/models/' + data0 + traindata + '_' + file_suffix + '.pkl')
explainer = shap.Explainer(loaded_model.predict, X_test)
shap_values = explainer(X_test)

from matplotlib import pyplot as plt
shap.summary_plot(shap_values, X_test, max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/traditional_model/summary/shap_' + data0 + traindata + 'p1.png', bbox_inches='tight')
plt.close()

shap.summary_plot(shap_values, X_test, plot_type="bar", max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/traditional_model/summary/shap_' + data0 + traindata + 'p2.png', bbox_inches='tight')
plt.close()


# Rename features in SHAP plots
mapping = pd.read_csv('/home/u7080475/ketamine/map_new.csv')
temp = shap_values

for i in range(len(temp.feature_names)):
    try:
        code_name = temp.feature_names[i]
        tempindx = list(mapping['CodeName2']).index(code_name)
        temp.feature_names[i] = mapping['Phenotype'][tempindx]
    except:
        pass

shap.summary_plot(temp, X_test, max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/traditional_model/summary/shap_'+data0+traindata+'p1_rename.png', dpi=600, bbox_inches='tight')
plt.close()

shap.summary_plot(temp, X_test, plot_type="bar", max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/traditional_model/summary/shap_'+data0+traindata+'p2_rename.png', bbox_inches='tight', dpi=600)
plt.close()




