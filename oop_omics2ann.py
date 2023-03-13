import argparse
import pandas as pd
import anndata as ad
import sys

class AnndataCreator:

    def __init__(self, mat_path, samp_path, umap_path, csv_file_path):
        self.mat_path = mat_path
        self.samp_path = samp_path
        self.umap_path = umap_path
        self.csv_file_path = csv_file_path
        self.mat_df = None
        self.samp_df = None
        self.umap_df = None
        self.type_dict = None
        self.adata = None

    def load_dataframes(self):
        self.mat_df = pd.read_csv(self.mat_path, index_col=0)
        self.samp_df = pd.read_csv(self.samp_path, index_col=0)
        self.umap_df = pd.read_csv(self.umap_path, index_col=0)

    def load_csv_file(self):
        self.type_dict = {}
        with open(self.csv_file_path, 'r') as f:
            for line in f:
                values = line.strip().split(',')
                field_name = values[0]
                string_type = values[1]
                self.type_dict[field_name] = string_type

    def check_field_names(self):
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
        for field in self.type_dict:
            self.samp_df[field] = self.samp_df[field].astype(self.type_dict[field])

    def transpose_matrix(self):
        mat_rows, mat_cols = self.mat_df.shape
        samp_rows, samp_cols = self.samp_df.shape
        if mat_cols == samp_rows:
            self.mat_df = self.mat_df.T

    def quality_check(self):
    # Check that the number of rows in the data matrix (X) matches the number of rows in the observation annotations (obs).
    # Check that the number of columns in the data matrix (X) matches the number of rows in the UMAP data frame (obsm).
    # Check that the row labels (index) in the data matrix (X) match the row labels in the observation annotations (obs).
    # Check that the row labels (index) in the data matrix (X) match the row labels in the UMAP data frame (obsm).
        if self.mat_df.shape[0] != self.samp_df.shape[0]:
            print("\nThe number of rows in X are:{}, The number of rows in obs is: {}".format(self.mat_df.shape[0], self.samp_df.shape[0]))
            raise ValueError("The number of rows in X and obs do not match.")
        if self.mat_df.shape[0] != self.umap_df.shape[0]:
            print("\nThe number of columns in X is:{}, The number of rows in obs is:{}".format(self.mat_df.shape[0], self.samp_df.shape[0]))
            raise ValueError("The number of columns in X and the number of rows in obsm do not match.")
        if not self.mat_df.index.equals(self.samp_df.index):
            print("\nThe number of rows in X are: {}, The number of rows in obs are: {}".format(self.mat_df.index, self.samp_df.index))
            raise ValueError("The row labels in X and obs do not match.\n")
        if not self.mat_df.index.equals(self.umap_df.index):
            print("\nThe number of rows in X are: {}, The number of rows in obsm are: {}".format(self.mat_df.index, self.samp_df.index))
            raise ValueError("The row labels in X and obsm do not match.\n")

    def create_anndata(self):
        self.adata = ad.AnnData(X=self.mat_df, obs=self.samp_df, dtype=self.mat_df.values.dtype)

    def check_umap_names(self):
        assert (self.umap_df.index == self.samp_df.index).all(), 'UMAP and sample metadata do not match'
        assert (self.umap_df.index == self.mat_df.index).all(), 'UMAP and count matrix do not match'

    def add_umap(self):
        umap_array = self.umap_df.to_numpy()
        self.adata.obsm['X_umap'] = umap_array 

    def write_h5ad_file(self, output_file_path):
        try:
            self.adata.write_h5ad(output_file_path)
        except Exception as e:
            print("\nCannot Write to h5ad because:{} ", e)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Process CSV files for AnnData object creation')
    parser.add_argument('-m', '--mat', required=True, help='File path for mat.csv')
    parser.add_argument('-s', '--samp', required=True, help='File path for samp.dat.csv')
    parser.add_argument('-u', '--umap', required=False, help='File path for umap.coord.csv')
    parser.add_argument('-c', '--csv', required=False, help='File path for example.csv')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path for h5ad file')
    args = parser.parse_args()

    # Create an AnndataCreator object
    creator = AnndataCreator(args.mat, args.samp, args.umap, args.csv)

    # Load dataframes
    creator.load_dataframes()

    # Load csv file
    creator.load_csv_file()

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

    # Check UMAP names
    creator.check_umap_names()

    # Add UMAP
    creator.add_umap()

    # Write h5ad file
    creator.write_h5ad_file(args.output)

if __name__ == "__main__":
    main()
