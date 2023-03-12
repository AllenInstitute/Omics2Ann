import argparse
import pandas as pd
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
import sys

# Define command line arguments
parser = argparse.ArgumentParser(description='Process CSV files for AnnData object creation')
parser.add_argument('--mat', type=str, required=True, help='File path for mat.csv')
parser.add_argument('--samp', type=str, required=True, help='File path for samp.dat.csv')
parser.add_argument('--umap', type=str, required=True, help='File path for umap.coord.csv')
args = parser.parse_args()

# Load the CSV files into pandas dataframes
mat_df = pd.read_csv(args.mat, index_col=0)
samp_df = pd.read_csv(args.samp, index_col=0)
umap_df = pd.read_csv(args.umap, index_col=0)

# Create an empty dictionary to store the metadata
type_dict = dict()

# Parse the example.csv file and store the metadata in type_dict
# Please modify line 20 to the name of the csv file
with open('example.csv', 'r') as f:
    # Iterate over the lines of the file
    for line in f:
        # Split each line into a list of values
        values = line.strip().split(',')
        # Use the first value as the key and the second value as the value
        field_name = values[0]
        string_type = values[1]
        # Add the key-value pair to the dictionary
        type_dict[field_name] = string_type

# Check if the fields in the metadata match the fields in the sample dataframe
samp_col_list = samp_df.columns.tolist()
type_dict_list = list(type_dict.keys())
diff_list = list(set(type_dict_list) - set(samp_col_list))
if len(diff_list) == 0:
    pass
else:
    print('Please check field in csv the file matches the field names in the metadata.')
    print('List of mismatching field names: {}'.format(diff_list))
    sys.exit()

# Convert the columns in the sample dataframe to the appropriate data types
for field in type_dict:
    samp_df[field] = samp_df[field].astype(type_dict[field])

# Transpose the count matrix if necessary
mat_rows, mat_cols = mat_df.shape
samp_rows, samp_cols = samp_df.shape
if mat_cols == samp_rows:
    mat_df = mat_df.T

# Define a function to perform quality checks on the data
def qualityCheck(X, obs, obsm):
    """
    Check that the number of rows in the data matrix (X) matches the number of rows in the observation annotations (obs).
    Check that the number of columns in the data matrix (X) matches the number of rows in the UMAP data frame (obsm).
    Check that the row labels (index) in the data matrix (X) match the row labels in the observation annotations (obs).
    Check that the row labels (index) in the data matrix (X) match the row labels in the UMAP data frame (obsm).
    """
    # Check if the number of rows in X matches the number of rows in obs
    if X.shape[0] != obs.shape[0]:
        print("\nThe number of rows in X are:{}, The number of rows in obs is: {}".format(X.shape[0], obs.shape[0]))
        raise ValueError("The number of rows in X and obs do not match.")

    # check if the number of columns in X matches the number of rows in obsm
    if X.shape[0] != obsm.shape[0]:
        print("\nThe number of columns in X is:{}, The number of rows in obs is:{}".format(X.shape[0], obs.shape[0]))
        raise ValueError("The number of columns in X and the number of rows in obsm do not match.")

    # check if the row labels in X match the row labels in obs
    if not X.index.equals(obs.index):
        print("\nThe number of rows in X are: {}, The number of rows in obs are: {}".format(X.index, obs.index))
        raise ValueError("The row labels in X and obs do not match.\n")

    # check if the row labels in X match the row labels in obsm
    if not X.index.equals(obsm.index):
        print("\nThe number of rows in X are: {}, The number of rows in obsm are: {}".format(X.index, obs.index))
        raise ValueError("The row labels in X and obsm do not match.\n")

def create_anndata(mat_th_df, samp_df):
    '''
    Create AnnData object to store expression and metadata
    '''
    adata = ad.AnnData(X=mat_th_df, obs=samp_df, dtype=mat_th_df.values.dtype)
    return adata

# Perform Quality Checks on the Data
qualityCheck(mat_df, samp_df, umap_df)

# Create the anndata object to store the expression and metadata
adata = create_anndata(mat_df, samp_df)

# Quality Checks for UMAP in obsm: Check that row names match
assert (umap_df.index == samp_df.index).all(), 'UMAP and sample metadata do not match'
assert (umap_df.index == mat_df.index).all(), 'UMAP and count matrix do not match'

umap_array = umap_df.to_numpy()
adata.obsm['X_umap'] = umap_array 

try:
    adata.write_h5ad('Mouse_medulla_10x_rna_rsc_319_integration.h5ad')
except Exception as e:
    print("\nCannot Write to h5ad because:{} ", e)
