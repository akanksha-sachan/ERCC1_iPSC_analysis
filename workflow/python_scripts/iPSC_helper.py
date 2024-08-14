import os
import sys
import numpy as np
import pandas as pd

import multiprocessing as mp

###### helper functions for the analysis of the iPSC data ######

###### Protein data processing ######

# matrix creation function for the protein data containing covariate samples in one matrix
def get_sample_by_protein_matrix(df, sample_identifier):
    """
    Filters the DataFrame for the given sample identifier, sorts by 'Sample Name', 
    selects specific columns, and keeps only the first gene name.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the sample data.
    sample_identifier (str): The string identifier for the sample (e.g., 'CM_KO_W2').
    
    Returns:
    pd.DataFrame: The processed DataFrame.
    """
    # Filter the DataFrame for the specific sample identifier
    filtered_df = df[df['Sample Name'].str.contains(sample_identifier)]
    # Sort the DataFrame by 'Sample Name'
    sorted_df = filtered_df.sort_values(by='Sample Name')
    # Select the relevant columns
    selected_df = sorted_df[['Sample Name', 'Gene Names', 'Intensity (Log10)', 'Normalized Intensity (Log10)']]
    # Keep only the first gene name if there are multiple gene names
    selected_df['Gene Names'] = selected_df['Gene Names'].str.split(';').str[0]
    # Print the first few rows and the shape of the processed DataFrame
    print(selected_df.head())
    print(selected_df.shape)
    
    return selected_df

def multiprocess_hdf5(df, sample_identifiers):
    pass
