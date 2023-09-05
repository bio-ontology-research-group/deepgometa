#!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from collections import Counter
from utils import Ontology, FUNC_DICT, NAMESPACES
import logging
import torch

logging.basicConfig(level=logging.INFO)

@ck.command()
@ck.option(
    '--data-root', '-dr', default='data',
    help='Prediction model')
@ck.option(
    '--go-file', '-gf', default='data/go-basic.obo',
    help='Gene Ontology file in OBO Format')
@ck.option(
    '--old-data-file', '-odf', default='data/swissprot_exp_2023_01.pkl',
    help='Uniprot KB, generated with uni2pandas.py')
@ck.option(
    '--new-data-file', '-ndf', default='data/swissprot_exp_2023_03.pkl',
    help='Uniprot KB, generated with uni2pandas.py')
def main(data_root, go_file, old_data_file, new_data_file):
    go = Ontology(go_file, with_rels=True)
    logging.info('GO loaded')
    
    df = pd.read_pickle(old_data_file)
    new_df = pd.read_pickle(new_data_file)
    
    print("OLD DATA FILE", len(df))
    print("NEW DATA FILE", len(new_df))
    
    logging.info('Processing annotations')

    annotations = list()

    for ont in ['cc', 'bp', 'mf']:
        out_terms_file = f'{data_root}/{ont}/terms.pkl'
        out_interpros_file = f'{data_root}/{ont}/interpros.pkl'
        train_data_file = f'{data_root}/{ont}/train_data.pkl'
        valid_data_file = f'{data_root}/{ont}/valid_data.pkl'
        test_data_file = f'{data_root}/{ont}/test_data.pkl'
        cnt = Counter()
        iprs = Counter()
        index = []
        for i, row in enumerate(df.itertuples()):
            ok = False
            for term in row.prop_annotations:
                if term != FUNC_DICT[ont] and go.get_namespace(term) == NAMESPACES[ont]:
                    cnt[term] += 1
                    ok = True
            if ok:
                for ipr in row.interpros:
                    iprs[ipr] += 1
                index.append(i)
            
        train_df = df.iloc[index]
        train_prots = set(train_df['proteins'])
        
        # Sort GO classes by their number of annotations
        terms = [g_id for g_id, c in sorted(cnt.most_common(), key=lambda x: (x[1], x[0]))]
        interpros = list(iprs.keys())

        print(f'Number of {ont} terms {len(terms)}')
        print(f'Number of {ont} iprs {len(iprs)}')

        terms_df = pd.DataFrame({'gos': terms})
        terms_df.to_pickle(out_terms_file)
        iprs_df = pd.DataFrame({'interpros': interpros})
        iprs_df.to_pickle(f'data/{ont}/interpros.pkl')
        
        # Split train/valid
        index = np.arange(len(train_df))
        np.random.shuffle(index)
        train_n = int(len(train_df) * 0.95)
        valid_df = train_df.iloc[index[train_n:]]
        valid_df.to_pickle(valid_data_file)
        train_df = train_df.iloc[index[:train_n]]
        train_df.to_pickle(train_data_file)

        # Test proteins
        index = []
        for i, row in enumerate(new_df.itertuples()):
            if row.proteins in train_prots:
                continue
            ok = False
            for term in row.prop_annotations:
                if term == 'GO:0005515' or term == FUNC_DICT[ont]:
                    continue
                if go.get_namespace(term) == NAMESPACES[ont]:
                    ok = True
                    break
            if ok:
                index.append(i)
        test_df = new_df.iloc[index]
        test_df.to_pickle(test_data_file)
        
        print(f'Train/Valid/Test proteins for {ont} {len(train_df)}/{len(valid_df)}/{len(test_df)}')


if __name__ == '__main__':
    main()
