import tempfile
from enum import Enum
from os.path import isfile, join
import dash_bootstrap_components as dbc
import base64
import os
import dash
import fire
import pandas as pd
import numpy as np
import logging

from appdata import AppDataPaths
from dash import Dash, dcc, html, no_update, Output, Input, MATCH, ALL
import plotly.graph_objects as go
from molcomplib import MolCompass
from molcompview.actions import init_callbacks, ColumnType
from molcompview.components import  molcompass_layout, header_alt


#Create enums for numerical, categorical and binary columns

def show(file,precompute=False,log_level="ERROR"):
    def make_dropdown(useful_columns):
        return dcc.Dropdown(
            id='dropdown',
            options=[{'label': i, 'value': i} for i in useful_columns],
            value=useful_columns[0]
        ) if len(useful_columns) > 0 else html.Div("No columns for coloring found")

    if not isfile(file):
        raise FileNotFoundError("File {} not found, please specify a csv file".format(file))

    logging.basicConfig(level=log_level)
    data = pd.read_csv(file)

    #Drop rows with NaN values in smiles x_coord or y_coord
    # data.dropna(subset=['smiles','x_coord','y_coord'],inplace=True)
    get_column('smiles',data,raise_error=True,modify_inplace=True)

    #Check, if there is a column with smiles
    if 'smiles' not in data.columns:
        raise ValueError("Dataframe should contain a column with smiles! Guess was unsuccessful. Please specify it with --smilesColumn")

    prob_column = get_column_prob(data)
    if prob_column != None:
        data.rename(columns={prob_column: 'logits'}, inplace=True)
        # get_column('class',data,raise_error=False,modify_inplace=True)
    try:
        x_col = get_column('x_coord',data,raise_error=True,modify_inplace=True)
        y_col = get_column('y_coord',data,raise_error=True,modify_inplace=True)
    except:
        logging.info("No x/y columns found, will calculate them by molcomplib")
        data = run_molcomplib(data,"smiles")
        data.rename(columns={'x': 'x_coord', 'y': 'y_coord'}, inplace=True)
        #Save the data with coordinates in a temporary directory
        temp_file = join(appath,file+"_temp.csv")
        data.to_csv(temp_file,index=False)


    #If we have both class and logits columns, we can calculate loss
    if 'class' not in data.columns:
        logging.error("No ground truth column found, MolCompass is running in limited mode. AD domain analysis will not be available")
    if 'class' in data.columns and 'logits' in data.columns:
        logging.error("We have both class and logits columns, calculating loss")
        data['loss'] = -data['class']*np.log(data['logits'])-(1-data['class'])*np.log(1-data['logits'])
    elif 'logits' in data.columns:
        logging.error("We have logits column, but no class column, calculating loss is not possible, AD domain analysis will not be available")
    else:
        logging.error("No logits column found, calculating loss is not possible, AD domain analysis will not be available")

    #Rename the columns x and y to x_coord and y_coord

    string_columns, categorical_columns, binary_columns, numerical_columns = _get_column_types(data)
    logging.info("String columns: {}".format(string_columns))
    logging.info("Categorical columns: {}".format(categorical_columns))
    logging.info("Binary columns: {}".format(binary_columns))
    logging.info("Numerical columns: {}".format(numerical_columns))
    #Convert to dictionary, where key is column name and value is the type of column
    column_types = {c: ColumnType.CATEGORICAL for c in data.columns if c in categorical_columns}
    column_types.update({c: ColumnType.BINARY for c in data.columns if c in binary_columns})
    column_types.update({c: ColumnType.NUMERICAL for c in data.columns if c in numerical_columns})
    useful_columns = categorical_columns + numerical_columns + binary_columns
    useful_columns = [col for col in useful_columns if col not in ["x_coord","X_coord","y_coord","Y_coord","smiles"]]
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.layout = dbc.Container([
        dbc.Row([dbc.Col(header_alt(), className='g-0')]),  # Header
        dbc.Row([dbc.Col(molcompass_layout(useful_columns), className='g-0')],id='main-layout'),  # Main
    ], fluid=True
    )
    init_callbacks(app,data,column_types)
    app.run_server(debug=True)


#
appath = None
def main():
    global appath
    app_paths = AppDataPaths()
    app_paths.setup()
    appath = app_paths.app_data_path
    print("Here appath is: ", appath)
    fire.Fire(show)

#if (__name__ == '__main__'):
#      print("Starting MolCompass")
      # fire.Fire(main)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
