import pandas as pd
import numpy as np

sheets = ['DeepDILI training set','DeepDILI test set','External test - NCTR','External test -Greene','External test - Xu','DrugBank']
models = ['mold2','mol2vec','maccs']
for model in models:
    for sheet in sheets[1:-1]:
        df = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheet)
        cols = ['Canonical SMILES','DILI_label',f'{model}_prob_DeepDILI']
        if sheet == 'DeepDILI test set':
            cols.append('initial_approval_year')
        df = df[cols]
        df = df.rename(columns={'Canonical SMILES': 'smiles', 'DILI_label': 'class'})
        df.rename(columns={f'{model}_prob_DeepDILI': 'prob'}, inplace=True)
        df.to_csv('../data/'+model+'_'+sheet+'.csv',index=False)

# data = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheets[1])
# data.to_csv("../data/deepdili_test.csv",index=False)
# data = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheets[2])
# data.to_csv("../data/deepdili_nctr.csv",index=False)
# data = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheets[3])
# data.to_csv("../data/deepdili_greene.csv",index=False)
# data = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheets[4])
# data.to_csv("../data/deepdili_xu.csv",index=False)
# data = pd.read_excel("../data/deepdili.xlsx",sheet_name=sheets[5])
# data.to_csv("../data/deepdili_drugbank.csv",index=False)
#
