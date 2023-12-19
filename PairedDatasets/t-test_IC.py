import pandas as pd
from scipy import stats
import sys
import os
import glob

def read_ic_values(file_path):
    """
    Reads the IC values from the file and returns them as a list of lists.
    Each sublist contains the IC values for one gene.
    """
    with open(file_path, 'r') as file:
        return [list(map(float, line.strip().split('\t'))) for line in file]

def perform_t_test(data):
    """
    Performs a t-test between the IC values of the first two genes in the data.
    Assumes that there are at least two genes in the data.
    """
    gene1_ic_values = data[0]
    gene2_ic_values = data[1]
    t_stat, p_value = stats.ttest_ind(gene1_ic_values, gene2_ic_values, equal_var=False)
    return p_value

# Path to the file with IC values
dir_path = ''  # Update this with your actual file path
files = glob.glob(os.path.join(dir_path, '*_IC.out'))

# Read the IC values
results = []
for file in files:
    ic_values = read_ic_values(file)
    p_value = perform_t_test(ic_values)
    results.append([os.path.basename(file), p_value])

df = pd.DataFrame(results, columns = ['file', 'p-value'])

df.to_csv('t-test_IC.tsv', sep = '\t', index = False)
