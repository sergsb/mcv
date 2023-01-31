from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.model_selection import cross_val_predict as cvp

_OPTIMAL_AMES_= {
 'max_depth': 8,
 'learning_rate': 0.1,
 'min_child_weight': 1,
 'gamma': 0.5,
 'subsample': 0.8,
 'colsample_bytree': 1.0}

_OPTIMAL_ER_= {
                 'max_depth': 8,
                 'learning_rate': 0.1,
                 'min_child_weight': 1,
                 'gamma': 0.5,
                 'subsample': 0.8,
                 'colsample_bytree': 1.0}

import rdkit
import pandas as pd
from rdkit import Chem
import numpy as np
import xgboost as xgb
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import molvs
from rdkit.Chem.SaltRemover import SaltRemover

def ecfp( mol, r=3, nBits=2048, errors_as_zeros=True):
    mol = Chem.MolFromSmiles(mol) if not isinstance(mol, rdkit.Chem.rdchem.Mol) else mol
    try:
        arr = np.zeros((1,))
        ConvertToNumpyArray(GetMorganFingerprintAsBitVect(mol, r, nBits), arr)
        return arr.astype(np.float32)
    except:
        return np.NaN if not errors_as_zeros else np.zeros((nBits,), dtype=np.float32)

ames = pd.read_csv("../data/Mutagenicity_N6512.csv")[['Canonical_Smiles','Activity']].rename(columns={'Canonical_Smiles':'smiles',"Activity":"y"})
er = pd.read_csv("../data/ER.csv")[['Smiles','Class']].rename(columns={'Smiles':'smiles',"Class":"y"})
s = molvs.Standardizer()
salts = SaltRemover()

def process(smiles): #Some of molecules are broken (2 instead of 1) we will remove this data
    try:
        m =  Chem.MolFromSmiles(smiles)
        m = s.standardize(m)
        m = salts(m)
        arr = np.zeros((1,))
        ConvertToNumpyArray(GetMorganFingerprintAsBitVect(m, 3, 2048), arr)
        return arr.astype(np.float32)
    except:
        return None


er['X'] = er['smiles'].apply(process)
ames['X'] = ames['smiles'].apply(process)
er.dropna(inplace=True)
ames.dropna(inplace=True)

#Build Xgboost CV models
def build_model(data,type):
    if type == 'ames':
        params = _OPTIMAL_AMES_
    else:
        params = _OPTIMAL_ER_
    # X_train, y_train, X_test, y_test = train_test_split(data,data,test_size=0.2,random_state=42)


    #Convert to DMatrix
    #data = xgb.DMatrix(np.vstack(data['X']), label=data['y'])
    #Define parameters
    params = {
        'n_estimators': 500,
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'tree_method': 'gpu_hist',
        'max_depth': params['max_depth'],
        'learning_rate': params['learning_rate'],
        'min_child_weight': params['min_child_weight'],
        'gamma': params['gamma'],
        'subsample': params['subsample'],
        'colsample_bytree': params['colsample_bytree'],
        'seed': 42
    }
    #Train model
    xgb_model = xgb.XGBClassifier(**params)
    y_pred = cvp(xgb_model, np.vstack(data['X']), data['y'], cv=5, n_jobs=1,verbose=10,method='predict_proba')
    #Cross validation
    #Train test split
    # X_train, X_test, y_train, y_test = train_test_split(data['X'], data['y'], test_size=0.2, random_state=42)
    # #Convert to DMatrix
    # dtrain = xgb.DMatrix(np.vstack(X_train['X']), label=y_train['y'])
    # dtest = xgb.DMatrix(np.vstack(X_test['X']), label=y_test['y'])
    # model = xgb.train(params, data, num_boost_round=1000, evals=[(dtest, 'test')], early_stopping_rounds=100, verbose_eval=10)
    # #Predict
    # y_pred = model.predict(dtest)
    # #Return dataframe of Smiles, y_true, y_pred
    # return pd.DataFrame({'smiles':X_test['smiles'],'y':y_test,'probs':y_pred})
    #Cross validation
    # cv_results = xgb.cv(
    #     params,
    #     data,
    #     num_boost_round=5000,
    #     nfold=5,
    #     metrics='auc',
    #     prediction=True,
    #     early_stopping_rounds=100,
    #     verbose_eval=10,
    #     seed=42)
    # return cv_results.pred
    #Save the predictions
    return pd.DataFrame({'smiles':data['smiles'].tolist(),'y':data['y'].tolist(),'probs':y_pred[:,1].tolist()})
build_model(ames,'er').to_csv('../data/er_predictions.csv',index=False)
    #np.save(f'../data/{type}_cv_predictions.npy', cv_results.pred)




