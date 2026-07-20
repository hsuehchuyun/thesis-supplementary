data0='FullW4_piHRSD_selected'

import pandas as pd
from sklearn.model_selection import train_test_split
ds = pd.read_csv('/home/u7080475/SSRI/'+data0+'.csv')
X= ds.drop(['ID','Site','piHRSDW4','as.factor(medication)fluoxetine','as.factor(medication)citalopram','as.factor(medication)paroxetine','as.factor(medication)escitalopram'], axis=1, errors='ignore')
y= ds['piHRSDW4']


from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest, f_regression

from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn import model_selection
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold
cv = model_selection.RepeatedKFold(n_splits=5, n_repeats=5, random_state=0) 

# ElasticNet 
from sklearn.linear_model import ElasticNet # ElasticNet
reg_model = ElasticNet()
param_grid = {'l1_ratio': [i * 0.1 for i in range(6,11)],'alpha': [10 ** i for i in range(-5, 3)], 'random_state': [0]}
reg_cv = GridSearchCV(reg_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40)
reg_cv.fit(X, y)
mse_folds = [
    reg_cv.cv_results_['split0_test_score'][reg_cv.best_index_],
    reg_cv.cv_results_['split1_test_score'][reg_cv.best_index_],
    reg_cv.cv_results_['split2_test_score'][reg_cv.best_index_],
    reg_cv.cv_results_['split3_test_score'][reg_cv.best_index_],
    reg_cv.cv_results_['split4_test_score'][reg_cv.best_index_]
]
print("LogReg each fold MSE:", mse_folds)
print("Best Hyperparameters:", reg_cv.best_params_)
print("Best Score:", -reg_cv.best_score_)


#KNN

from sklearn.neighbors import KNeighborsRegressor
knn = KNeighborsRegressor()
param_grid = {'n_neighbors': [3, 5, 10, 15, 20, 25, 30, 40, 50, 60], 'weights': ['uniform', 'distance'], 'p': [1, 2]}
knn_cv = GridSearchCV(knn, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40)
knn_cv.fit(X, y)
knn_mse_folds = [
    knn_cv.cv_results_['split0_test_score'][knn_cv.best_index_],
    knn_cv.cv_results_['split1_test_score'][knn_cv.best_index_],
    knn_cv.cv_results_['split2_test_score'][knn_cv.best_index_],
    knn_cv.cv_results_['split3_test_score'][knn_cv.best_index_],
    knn_cv.cv_results_['split4_test_score'][knn_cv.best_index_]
]
print("KNN each fold MSE:", knn_mse_folds)
print("Best Hyperparameters:", knn_cv.best_params_)
print("Best Score:", -knn_cv.best_score_)


#SVR
from sklearn.svm import SVR
svm = SVR()
param_grid = {'C': [10 ** i for i in range(2, -6, -1)], 'kernel': ['linear'], 'kernel': ['linear', 'rbf'], 'epsilon': [0.01, 0.1, 0.2, 0.5],'gamma': ['scale', 'auto']}
svm_cv = GridSearchCV(svm, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40)

svm_cv.fit(X, y)
svm_mse_folds = [
    svm_cv.cv_results_['split0_test_score'][svm_cv.best_index_],
    svm_cv.cv_results_['split1_test_score'][svm_cv.best_index_],
    svm_cv.cv_results_['split2_test_score'][svm_cv.best_index_],
    svm_cv.cv_results_['split3_test_score'][svm_cv.best_index_],
    svm_cv.cv_results_['split4_test_score'][svm_cv.best_index_]
]
print("SVM each fold MSE:", svm_mse_folds)
print("Best Hyperparameters:", svm_cv.best_params_)
print("Best Score:", -svm_cv.best_score_)
svmScore= -svm_cv.best_score_


#xgboost https://blog.csdn.net/wzmsltw/article/details/50994481
from xgboost import XGBRegressor
xgboost_model = XGBRegressor(nthread=1) 

param_grid = {
     	'n_estimators': [100], 
       	'booster': ['gbtree'], #'gblinear' 
        'learning_rate': [0.1], 
        'gamma': [0], 
        'reg_alpha': [0], 
        'reg_lambda': [0], 
        'objective':['reg:squarederror'],
        'colsample_bylevel': [0.7,0.9], 
        'colsample_bynode': [0.1,0.3], 
        'colsample_bytree': [0.1,0.3], 
        'subsample': [0.8],
        'min_child_weight': [5],  
        'max_depth': [1], 
        'random_state': [0]
        }
xgboost_cv = GridSearchCV(xgboost_model, param_grid, cv=cv, scoring='neg_mean_squared_error', n_jobs=40)  

xgboost_cv.fit(X, y) 
xgb_mse_folds = [
    xgboost_cv.cv_results_['split0_test_score'][xgboost_cv.best_index_],
    xgboost_cv.cv_results_['split1_test_score'][xgboost_cv.best_index_],
    xgboost_cv.cv_results_['split2_test_score'][xgboost_cv.best_index_],
    xgboost_cv.cv_results_['split3_test_score'][xgboost_cv.best_index_],
    xgboost_cv.cv_results_['split4_test_score'][xgboost_cv.best_index_]
]
print("XGBoost each fold mse:", xgb_mse_folds)
print("Best Hyperparameters:", xgboost_cv.best_params_)
print("Best Score:", -xgboost_cv.best_score_)


#ensemble
from sklearn.ensemble import VotingRegressor
ensem = [('elasticnet', reg_cv.best_estimator_),('knn', knn_cv.best_estimator_),('svm', svm_cv.best_estimator_), ('xgboost', xgboost_cv.best_estimator_)]
ensemble = VotingRegressor(ensem)
ensemble.fit(X, y)
print("Test R2 Score:", ensemble.score(X, y))

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score

res = pd.DataFrame([
    ['ElasticNet', -reg_cv.best_score_,
     mean_squared_error(y, reg_cv.predict(X)),
     r2_score(y, reg_cv.predict(X))],
    
    ['KNN', -knn_cv.best_score_,
     mean_squared_error(y, knn_cv.predict(X)),
     r2_score(y, knn_cv.predict(X))],
    
    ['SVR', -svm_cv.best_score_,
     mean_squared_error(y, svm_cv.predict(X)),
     r2_score(y, svm_cv.predict(X))],
    
    ['XGBoost', -xgboost_cv.best_score_,
     mean_squared_error(y, xgboost_cv.predict(X)),
     r2_score(y, xgboost_cv.predict(X))],
    
    ['Ensemble', 'NA',
     mean_squared_error(y, ensemble.predict(X)),
     r2_score(y, ensemble.predict(X))]
], columns=['Method', 'CV_MSE', 'MSE', 'R2']) 


print(res)
res.to_csv('/home/u7080475/SSRI/ML_cv/results/res'+data0+'.csv', index=False)


import pickle
save_dir = '/home/u7080475/SSRI/ML_cv/models/'
pickle.dump(reg_cv, open(save_dir + data0 +  '_elasticnet.pkl', 'wb'))
pickle.dump(knn_cv, open(save_dir + data0 +  '_knn.pkl', 'wb'))
pickle.dump(svm_cv, open(save_dir + data0 + '_svr.pkl', 'wb'))
pickle.dump(xgboost_cv, open(save_dir + data0 +  '_xgboost.pkl', 'wb'))
pickle.dump(ensemble, open(save_dir + data0 +  '_ensemble.pkl', 'wb'))


#========for SHAP
from multiprocessing import Pool
import pickle
import shap
import joblib

indx = res['MSE'].idxmin()
best_method = res['Method'][indx]
print("best method(SHAP):", best_method)

file_suffix = best_method.lower()
loaded_model = joblib.load('/home/u7080475/SSRI/ML_cv/models/' + data0  + '_' + file_suffix + '.pkl')
explainer = shap.Explainer(loaded_model.predict, X)
shap_values = explainer(X)

from matplotlib import pyplot as plt
shap.summary_plot(shap_values, X, max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/ML_cv/summary/shap_' + data0  + 'p1.png', bbox_inches='tight')
plt.close()

shap.summary_plot(shap_values, X, plot_type="bar", max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/ML_cv/summary/shap_' + data0 + 'p2.png', bbox_inches='tight')
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

shap.summary_plot(temp, X, max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/ML_cv/summary/shap_'+data0+'p1_rename.png', dpi=600, bbox_inches='tight')
plt.close()

shap.summary_plot(temp, X, plot_type="bar", max_display=20, show=False)
plt.savefig('/home/u7080475/SSRI/ML_cv/summary/shap_'+data0+'p2_rename.png', bbox_inches='tight', dpi=600)
plt.close()
