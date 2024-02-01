import pandas as pd
import os
import glob
import sys
import argparse

parser = argparse.ArgumentParser(description='Assign functions to genera')
parser.add_argument('--pkl_dir', required = True, help = 'Path to directory with Pkl files')
parser.add_argument('--accessions', required = True, help = 'Accessions for genera of interest')
parser.add_argument('--otu_table', required = True, help = "OTU table to be used for abundance")

args = parser.parse_args()

# Directory with the .pkl files
directory_path = args.pkl_dir

# Load the list of accession names from the txt file
accession_file_path = args.accessions
accession_df = pd.read_csv(accession_file_path, sep = '\t')
accession_names = accession_df['Accession'].to_list()

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

# Now unique_go_inds_dict contains all the unique go_ind values for each accession and type
unique_accessions = sorted({key.rsplit('_', 1)[0] for key in unique_go_inds_dict.keys()})
go_terms_df = pd.DataFrame(index=unique_accessions, columns=['cc', 'bp', 'mf'])
for key, go_list in unique_go_inds_dict.items():
    accession, file_type = key.rsplit('_', 1)
    go_terms_df.at[accession, file_type] = go_list
go_terms_df.to_csv('GO_terms.csv')

for index, row in  go_terms_df.iterrows():
    if row.isnull().any():  # Check if any value in the row is NaN
        print(f"Row {index} contains NaN values and will be deleted.")
        go_terms_df.drop(index, inplace=True)

##############

# Load mapping files
cc_mappings = pd.read_pickle(args.pkl_dir + 'terms_cc.pkl')
mf_mappings = pd.read_pickle(args.pkl_dir + 'terms_mf.pkl')
bp_mappings = pd.read_pickle(args.pkl_dir + 'terms_bp.pkl')

def map_indices_to_go(row, mapping_df):
    valid_indices = [int(i) for i in row if pd.notnull(i)]     # Filter out NaN values and convert to integers, but why are they there????
    return [mapping_df.iloc[i]['gos'] for i in valid_indices]

go_terms_df['cc'] = go_terms_df['cc'].apply(map_indices_to_go, mapping_df=cc_mappings)
go_terms_df['bp'] = go_terms_df['bp'].apply(map_indices_to_go, mapping_df=bp_mappings)
go_terms_df['mf'] = go_terms_df['mf'].apply(map_indices_to_go, mapping_df=mf_mappings)

def flatten_gos(row):
    return [gos for sublist in row for gos in sublist if gos is not None]

go_terms_df['gos'] = go_terms_df.apply(flatten_gos, axis=1)
go_terms_df = go_terms_df.drop(['cc', 'bp', 'mf'], axis=1)

all_go_terms = list(set(go_term for sublist in go_terms_df['gos'] for go_term in sublist)) # Flatten and join all GO terms

# Initialize a dictionary to store the GO term occurrences
go_occurrences = {taxa: {go_term: 0 for go_term in all_go_terms} for taxa in go_terms_df.index}
for taxa, go_list in go_terms_df['gos'].items():
    for go_term in go_list:
        go_occurrences[taxa][go_term] = 1
go_freq = pd.DataFrame.from_dict(go_occurrences, orient='index')

# Change GCA/GCF accessions to actual taxa names
mapping_dict = accession_df.set_index('Accession')['Genus'].to_dict()
go_freq.index = go_freq.index.map(mapping_dict)
go_freq.reset_index(inplace=True)


##############

# From OTU table to OTU dict
otu_table = pd.read_csv(args.otu_table, sep = '\t', header = 0)
otu_table.set_index('OTU', inplace=True)
otu_table_t = otu_table.T
otu_dict = otu_table_t.apply(lambda x: x.to_dict(), axis=1) # Convert each row to a dictionary with taxa as keys and abundances as values
otu_dict_df = pd.DataFrame(otu_dict, columns=['TaxaAbundance'])
otu_dict_df.reset_index(inplace=True)
otu_dict_df.rename(columns={'index': 'Sample'}, inplace=True)

##############

go_freq.set_index('index', inplace=True)
taxa_abundance_df = pd.DataFrame(list(otu_dict_df['TaxaAbundance'])).set_index(otu_dict_df['Sample'])
go_function_abundance_df = pd.DataFrame(index=all_go_terms, columns=otu_dict_df['Sample'])

for go_term in all_go_terms: # Iterate over each GO term to calculate its abundance in each sample
    for sample in otu_dict_df['Sample']: # For each sample, calculate the abundance of the GO term
        go_function_abundance_df.at[go_term, sample] = (go_freq[go_term] * taxa_abundance_df.loc[sample]).sum() # Multiply the taxa presence/absence by their abundance and sum


go_function_abundance_df.to_csv("OTU_function_abundance.tsv", sep = "\t")
