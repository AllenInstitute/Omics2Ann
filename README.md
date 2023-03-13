Omics2Ann
=================================================
![cover](images/sheet_to_text.png)

## About Omics2Ann
Omics data, such as single-cell RNA sequencing data, contains vast amounts of information that can provide insight into biological processes at the cellular level. However, analysis of this data requires the integration of multiple datasets, each containing different types of information, such as expression levels and metadata. The code presented here provides a framework for integrating multiple CSV files containing gene expression and metadata information into an AnnData object, which can be converted into a h5ad file. The h5ad file can then be uploaded onto Cirrocumulus, a web-based platform for single-cell RNA sequencing data analysis, for further processing and visualization.

Cirrocumulus provides a user-friendly interface for exploring the data and identifying patterns that may not be apparent from individual datasets. It also offers a range of tools for analyzing the data, including clustering, differential gene expression, and pathway analysis. By uploading the h5ad file generated by this code onto Cirrocumulus, users can take advantage of these features to gain a more comprehensive understanding of the underlying biological processes.

Overall, the integration of multiple datasets and the use of web-based platforms like Cirrocumulus can provide valuable insights into complex biological systems. This can lead to a deeper understanding of the data, identification of novel biological processes, and new research questions that can be pursued. By exploring the data in a comprehensive manner, researchers can gain insights that would have been missed when analyzing individual datasets. This, in turn, can facilitate the generation of new hypotheses and inform future experimental designs. Therefore, the use of this code and Cirrocumulus platform can contribute significantly to the advancement of our understanding of complex biological systems.

Table of Contents
-----------------
* [Usage](#usage)
* [Authors and history](#authors-and-history)
* [Acknowledgments](#acknowledgments)
* [References](#references)

## Start Guide
1. Download the files from Github
```bash
https://github.com/beagan-svg/Omics2Ann
```
2. Install the necessary Python packages using Conda
  - Open a terminal and copy packages.yml into a working directory
  - Navigate to the directory containing packages.yml
  - Run the following command to create a new conda environment named pyAnn
```
conda env create -f packages.yml -n pyAnn
```
3. Once the environment is created, activate it by running the following command
```
conda activate pyAnn
```
4. Copy create_anndata.py into the working directory and modify the line code containing and change to desired name
```
adata.write_h5ad('Mouse_medulla_10x_rna_rsc_319_integration.h5ad') 
```
5. Run the script by running the following command to generate the h5ad file
```
python3 create_anndata.py -m path/to/mat.csv -s path/to/samp.dat.csv -u path/to/umap.coord.csv
```
6. Convert h5ad to cirro format by running the following command. Make sure your in the same directory as the generated h5ad file
```
cirro prepare_data --format parquet --no-auto-groups file_name.h5ad
```
7. Follow this page [confluence](http://confluence.corp.alleninstitute.org/pages/viewpage.action?spaceKey=BIOIN&title=Cirrocumulus) to upload onto Cirrocumulus

## Authors and History

* Beagan Nguy - Algorithm Design
* Anish Chakka - Project Manager

## Acknowledgments

Allen Institute Bioinformatics Core Team

