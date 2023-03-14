import argparse
import pandas as pd
import anndata as ad
import sys
import os

class AnndataCreator:

    def __init__(self, mat_path, samp_path, umap_path = None, data_type_path = None):
        '''
        Initializes the paths for the expression, metadata, umap, and data type file
        '''
        self.mat_path = mat_path
        self.samp_path = samp_path
        self.umap_path = umap_path
        self.data_type_path = data_type_path
        self.mat_df = None
        self.samp_df = None
        self.umap_df = None
        self.type_dict = None
        self.adata = None

    def load_dataframes(self):
        '''
        Load inputs into pandas data frame
        '''
        self.mat_df = pd.read_csv(self.mat_path, index_col=0)
        self.samp_df = pd.read_csv(self.samp_path, index_col=0)
        self.umap_df = pd.read_csv(self.umap_path, index_col=0)

    def load_data_type_file(self):
        '''
        Load data type csv file into a dictionary
        '''
        self.type_dict = {}
        with open(self.data_type_path, 'r') as f:
            # Iterate over the lines of the file
            for line in f:
                # Split each line into a list of values
                values = line.strip().split(',')
                # Key = field name, Value = data type
                field_name = values[0]
                string_type = values[1]
                self.type_dict[field_name] = string_type

    def check_field_names(self):
        '''
        Checks if the columns in the metadata data frame matches the field names in data type csv file.
        '''
        samp_col_list = self.samp_df.columns.tolist()
        type_dict_list = list(self.type_dict.keys())
        diff_list = list(set(type_dict_list) - set(samp_col_list))
        if len(diff_list) == 0:
            pass
        else:
            print('Please check field in csv the file matches the field names in the metadata.')
            print('List of mismatching field names: {}'.format(diff_list))
            sys.exit()

    def cast_columns(self):
        '''
        Casts the columns in the metadata data frame to their corresponding data types as specified in the data types dictionary
        '''
        if self.type_dict != None:
            for field in self.type_dict:
                self.samp_df[field] = self.samp_df[field].astype(self.type_dict[field])

    def transpose_matrix(self):
        '''
        Transpose expression data matrix if sample and gene names are switch between rows and columns
        This is necessary to match X, obs when converting to AnnData
        '''
        mat_rows, mat_cols = self.mat_df.shape
        samp_rows, samp_cols = self.samp_df.shape
        if mat_cols == samp_rows:
            self.mat_df = self.mat_df.T

    def quality_check(self):
        # Check that the number of rows in the data matrix (X) matches the number of rows in the observation annotations (obs).
        # Check that the row labels (index) in the data matrix (X) match the row labels in the observation annotations (obs).
        if self.mat_df.shape[0] != self.samp_df.shape[0]:
            print("\nThe number of rows in X are:{}, The number of rows in obs is: {}".format(self.mat_df.shape[0], self.samp_df.shape[0]))
            raise ValueError("The number of rows in X and obs do not match.")
        if not self.mat_df.index.equals(self.samp_df.index):
            print("\nThe number of rows in X are: {}, The number of rows in obs are: {}".format(self.mat_df.index, self.samp_df.index))
            raise ValueError("The row labels in X and obs do not match.\n")

    def create_anndata(self):
        '''
        Creates an AnnData object with the loaded expression and metadata data frames. The created object is stored in the adata variable
        '''
        print(self.mat_df)
        print(self.samp_df)
        self.adata = ad.AnnData(X=self.mat_df, obs=self.samp_df, dtype=self.mat_df.values.dtype)

    def check_umap_names(self):
        '''
        Check that the row labels (index) in the UMAP dataframe match the row labels (index) in the expression and metadata dataframe.
        '''
        umap_index = self.umap_df.index
        samp_index = self.samp_df.index
        mat_index = self.mat_df.index
        assert umap_index.equals(samp_index), f"Row labels (index) in UMAP dataframe do not match row labels in metadata dataframe. UMAP index: {umap_index}, metadata dataframe: {samp_index}"
        assert umap_index.equals(mat_index), f"Row labels (index) in UMAP dataframe do not match row labels in count matrix. UMAP index: {umap_index}, expression count dataframe: {mat_index}"

    def add_umap(self):
        '''
        Adds the UMAP coordinates as a variable in the obsm attribute of the adata object.
        '''
        print(self.umap_df)
        umap_array = self.umap_df.to_numpy()
        self.adata.obsm['X_umap'] = umap_array 

    def write_h5ad_file(self, output_file_path):
        '''
        Writes the adata object to an h5ad file at the specified output file path. If the write operation fails, it prints an error message with the reason for the failure.
        '''
        try:
            self.adata.write_h5ad(output_file_path)
            print("h5ad file has been generated here: {}".format(os.getcwd()))
        except Exception as e:
            print("\nCannot Write to h5ad because:{} ", e)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Process expression count matrix, metadata, and umap files for AnnData object creation')
    usage = '''
        %(prog)s -m <expression_count_matrix_file_path> -s <metadata_file_path> -o <output_h5ad_file_path> [-u <umap_file_path>] [-c <data_types_file_path>]
    '''
    parser.usage = usage
    parser.add_argument('-m', '--mat', required=True, help='File path for mat.csv (expression counts matrix)')
    parser.add_argument('-s', '--samp', required=True, help='File path for samp.dat.csv (metadata)')
    parser.add_argument('-u', '--umap', required=False, help='File path for umap.coord.csv (umap coordinates)')
    parser.add_argument('-dt', '--dtypes', required=False, help='File path for data_types.csv')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path for h5ad file. Remember to include the .h5ad extension to the file name')
    parser.add_argument('-wd', '--setwd', type=str, required=False, help='Set the work directory')
    args = parser.parse_args()
    #print(vars(args))

    # Check if working directory exist
    if args.setwd:
        try:
            os.chdir(args.setwd)
        except:
            raise SystemExit("\n{} (working directory) does not exist. Please check your working directory pathway again.\n".format(args.setwd))

    # Create an AnndataCreator object
    creator = AnndataCreator(args.mat, args.samp, args.umap, args.dtypes)

    # Load dataframes
    creator.load_dataframes()

    # Load data type file
    creator.load_data_type_file()

    # Check field names
    creator.check_field_names()

    # Cast columns
    creator.cast_columns()

    # Transpose matrix
    creator.transpose_matrix()

    # Quality check
    creator.quality_check()

    # Create anndata
    creator.create_anndata()
    
    if args.umap:
        # Check UMAP names
        creator.check_umap_names()

        # Add UMAP
        creator.add_umap()

    # Write h5ad file
    creator.write_h5ad_file(args.output)

if __name__ == "__main__":
    main()
