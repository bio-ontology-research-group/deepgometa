import argparse
import pandas as pd

# Function to extract the accession from the path
def extract_accession(path):
    parts = path.strip().split("\\")
    return parts[-1][:15] if len(parts) > 1 else None

# Set up argument parser
parser = argparse.ArgumentParser(description='Extract accessions from OTU table and rdp-genomes_clean.txt')
parser.add_argument('--otu_table', required=True, help='Path to OTU table file')
parser.add_argument('--rdp_genomes_clean', required=True, help='Path to rdp-genomes_clean.txt file')
parser.add_argument('--output', required=True, help='Path to output file')

# Parse arguments
args = parser.parse_args()

# Read the OTU table
otu_df = pd.read_csv(args.otu_table, sep='\t')
genus_names = otu_df.iloc[:, 0].unique()  # Extract unique genus names
processed_genus_names = [name.split('_')[0] for name in genus_names]

# Read the rdp-genomes_clean.txt file and create a mapping
genus_to_accession = {}
with open(args.rdp_genomes_clean, 'r') as file:
    for line in file:
        genus, path = line.strip().split('\t')
        accession = extract_accession(path)
        genus_to_accession[genus] = accession

# Create a list for new DataFrame rows
new_rows = []
for genus in processed_genus_names:
    if genus in genus_to_accession:
        new_rows.append({'Genus': genus, 'Accession': genus_to_accession[genus]})
    else:
        print(f"Genus '{genus}' not found in rdp-genomes_clean.txt")

# Create a new DataFrame from the list of rows
new_df = pd.DataFrame(new_rows, columns=['Genus', 'Accession'])

# Save the DataFrame to the specified output file
new_df.to_csv(args.output, sep='\t', index=False)
