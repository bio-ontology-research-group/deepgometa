#!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
import gzip
import os
from collections import Counter
import logging

logging.basicConfig(level=logging.INFO)

@ck.command()
def main():
    # Load interpro data
    df1 = pd.read_pickle('data/swissprot_exp_2022_04.pkl')
    #df3 = pd.read_pickle('data/swissprot_exp_2023_03.pkl')
    seqs = {}
    for df in [df1,]:
        for row in df.itertuples():
            seqs[row.proteins] = row.sequences
    
    with open('data/esm_missing.fa', 'w') as f:
        for prot_id, seq in seqs.items():
            org = prot_id.split('_')[1]
            if os.path.exists(f'data/esm15B/{org}/{prot_id}.pt'):
                continue
            f.write('>' + prot_id + '\n')
            f.write(seq + '\n')

if __name__ == '__main__':
    main()
