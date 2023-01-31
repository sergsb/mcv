from ray.tune.sklearn import TuneGridSearchCV
from xgboost import XGBClassifier
from xgboost_ray import RayDMatrix, RayParams, train, RayXGBClassifier
import rdkit
import pandas as pd
from rdkit import Chem
import numpy as np
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


#%%
er['X'] = er['smiles'].apply(process)
ames['X'] = ames['smiles'].apply(process)
er.dropna(inplace=True)
ames.dropna(inplace=True)



ray_params = RayParams(
    num_actors=8,
    gpus_per_actor=1,
    cpus_per_actor=4,   # Divide evenly across actors per machine
)
def gridsearch(data):
    estimator = XGBClassifier(
        objective='binary:logistic',
        # resources_per_trial={'gpu': 1},
        tree_method='gpu_hist',
        eval_metric='auc',
        n_jobs=4,
        seed=42
    )
    parameters = {
        'n_estimators': [100],
        'max_depth': range(2, 10, 2),
        'learning_rate': [0.1, 0.01, 0.05],
        'min_child_weight': [1, 5, 10],
        'gamma': [0.5, 1, 1.5, 2, 5],
        'subsample': [0.6, 0.8, 1.0],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'ray_params': [ray_params]
    }
    grid_search = TuneGridSearchCV(
        estimator=estimator,
        param_grid=parameters,
        scoring='roc_auc',
        n_jobs=1,
        cv=5,  #predefied_split(),
    )
    return grid_search.fit(np.vstack(data['X'].values),np.vstack(data['y'].values))

res = gridsearch(er)