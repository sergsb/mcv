import hashlib
from os.path import isfile, join
import logging
import pandas as pd
from appdata import AppDataPaths

#Make loggers
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app_paths = AppDataPaths()
app_paths.setup()

def _get_column_types(data):
    approp_columns = list(data.select_dtypes(include=['object','float','double','int']).columns)
    logging.info("String columns: {}".format(approp_columns))
    #Select all categorical columns
    categorical_columns = [col for col in approp_columns if 2 < data[col].nunique() < 10]
    logging.info("Categorical columns: {}".format(categorical_columns))
    binary_columns = [col for col in approp_columns if data[col].nunique() == 2]
    logging.info("Binary columns: {}".format(binary_columns))
    #Select all non-categorical columns with numerical values, check if they are double
    numerical_columns = [col for col in data.columns if col in approp_columns and data[col].dtype == float]
    logging.info("Numerical columns: {}".format(numerical_columns))
    return approp_columns, categorical_columns, binary_columns, numerical_columns

def get_column_(data, list_of_column_names,name,raise_error=True,modify_inplace=False):
    column = [x for x in data.columns if x.lower() in list_of_column_names]
    if raise_error:
        assert len(column) == 1, "Dataframe should contain ONLY one {} column, but found: {}".format(name,column)
    else:
        if len(column) == 0:
            return None
    if modify_inplace:
        data.rename(columns={column[0]: name}, inplace=True)
    return column[0]

def get_column(column,data,raise_error=True,modify_inplace=False):
    correspondence = {'smiles':['smiles', 'smiles', 'smiles', 'molecules', 'structures', 'mols', 'smi','canonical_smiles','canonical_smi','canonicalsmiles','canonical smiles'],
                      'class':['class', 'classes', 'active','act','target','targets'],
                      'train/test split':['train','test','split','train/test','set'],
                      'logits':['logits','prob'],
                      'x_coord':['x_coord','X_coord'],
                      'y_coord':['y_coord','Y_coord']}
    return get_column_(data,correspondence[column],column,raise_error,modify_inplace)

def run_molcomplib(data,smilesColumn):
    compass = MolCompass()
    #Check if there is y column
    if 'y' in data.columns:
        data.rename(columns={'y': 'class'}, inplace=True)
    data = compass.process(data)
    return data


def get_column_prob(data):
    guess = get_column('logits',data,raise_error=False)
    if guess == None:
        logging.info("No logits column found with standard names, trying to guess the probability column")
        #Try to guess the column with probabilities
        #Get the distribution of values in all columns, and min and max
        min_max = data.describe().loc[['min','max']]
        #Get the columns with values between 0 and 1 but not 0 or 1
        prob_columns = [col for col in min_max.columns if (min_max[col]['min'] >= 0 and min_max[col]['max'] <= 1 and min_max[col]['min'] != 0 and min_max[col]['max'] != 1)]
        if len(prob_columns) == 1:
            guess = prob_columns[0]
            logging.info("Guessed column with probabilities: {}".format(guess))
        elif len(prob_columns) > 1:
            logging.info("More than one column with values between 0 and 1 found, please specify the column with probabilities")
            return None
        else:
            logging.info("Could not guess the column with probabilities, please specify it manually")
    return guess

def read_processed_file(filename):
    pass


def process_new_file(filename, smilesColumn, processed_file):
    data = pd.read_csv(filename, sep=',', header=0, index_col=False, low_memory=False)
    #Check if the file has a smiles column
    if smilesColumn == None:
        smilesColumn = get_column('smiles',data,raise_error=True)
        #Rename the column to smiles
        data.rename(columns={smilesColumn: 'smiles'}, inplace=True)

#    else:
#        data.rename(columns={smilesColumn: 'smiles'}, inplace=True)





def file_processing_entrypoint(filename,smilesColumn):
    #Check that the file exists and is a csv file
    if not isfile(filename):
        raise FileNotFoundError("File {} not found, please specify a csv file".format(filename))
    if not filename.endswith(".csv"):
        raise ValueError("File {} is not a csv file, please specify a csv file".format(filename))
    #Calculate md5 hash of the file
    hash = hashlib.md5(open(filename,'rb').read()).hexdigest()
    #Check if the file has been processed before, it is stored in app_paths.app_data_path
    processed_file = join(app_paths.app_data_path,hash+".csv")
    if not isfile(processed_file):
        process_new_file(filename,smilesColumn,processed_file)
    return read_processed_file(processed_file)

file_processing_entrypoint("../data/er_predictions.csv",None)


    # #Check if the first line in a file does not start from "Processed by MolCompass v 1"
    # with open(filename) as f:
    #     first_line = f.readline()
    #     if first_line.startswith("Processed by MolCompass v 1"):
    #         return read_processed_file(filename)
    #         #raise ValueError("File {} has already been processed by MolCompass, please specify a file that has not been processed yet".format(data))
