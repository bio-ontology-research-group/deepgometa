import pandas as pd
import os
import glob
import sys

dataset = ''

# Directory with the .pkl files
directory_path = ''

# Load the list of bacterial genome accession names from the txt file
accession_file_path = ''
with open(accession_file_path, 'r') as file:
    accession_names = [line.strip() for line in file]

unique_go_inds_dict = {}

for accession in accession_names:
    file_paths = glob.glob(f"{directory_path}/{accession}*_preds_*.pkl")
    if not file_paths:
        print(f"No files found for accession: {accession}")
        continue
    for file_path in file_paths:
        df = pd.read_pickle(file_path)
        go_inds_list = df['go_ind'].explode().unique().tolist()
        file_type = file_path.split('_')[-1].split('.')[0]        
        key = f"{accession}_{file_type}"
        unique_go_inds_dict[key] = go_inds_list

##############

# Load mapping files
cc_mappings = pd.read_pickle('terms_cc.pkl')
mf_mappings = pd.read_pickle('terms_mf.pkl')
bp_mappings = pd.read_pickle('terms_bp.pkl')

def get_dataframe(key):
    _, term_type = key.split('_')
    if term_type.lower() == 'cc':
        return cc_mappings
    elif term_type.lower() == 'mf':
        return mf_mappings
    elif term_type.lower() == 'bp':
        return bp_mappings
    else:
        raise ValueError("Unknown term type")

def map_terms(example_dict, cc_mappings, mf_mappings, bp_mappings):
    mapped_dict = {}
    for key, indices in example_dict.items():
        dataframe = get_dataframe(key)
        terms_list = []
        for index in indices:
            if not pd.isna(index):
                terms_list.append(dataframe.loc[index, 'gos'])
        mapped_dict[key] = terms_list  
    return mapped_dict

mapped_dict = map_terms(unique_go_inds_dict, cc_mappings, mf_mappings, bp_mappings)

flat_mapped_dict = {}
for accession in accession_names:
    flat_mapped_dict[accession] = []

for key, val in mapped_dict.items():
    sample_parts = key.split('_')
    sample = sample_parts[0]
    flat_mapped_dict[sample].append(val)

flat_mapped_dict = {key: [item for sublist in value for item in sublist] for key, value in flat_mapped_dict.items()}

##############

all_functions = []
for key, functions in flat_mapped_dict.items():
    for function in functions:
        all_functions.append(function)

all_functions = list(set(all_functions))

binary_df = pd.DataFrame(index = accession_names, columns = all_functions)
binary_df[:] = 0

for sample, go_functions in flat_mapped_dict.items():
    for go in go_functions:
        if go in binary_df.columns:
            binary_df.at[sample, go] = 1

binary_df.to_csv(directory_path + '/binary_function_df.tsv', sep = '\t')
