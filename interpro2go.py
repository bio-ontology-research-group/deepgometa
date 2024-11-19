#!/usr/bin/env python

import numpy as np
import pandas as pd
import click as ck
import sys
from collections import deque, Counter
import time
import logging
from utils import FUNC_DICT, Ontology, NAMESPACES

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@ck.command()
@ck.option(
    '--data-root', '-dr', default='data',
    help='Prediction model')
@ck.option(
    '--ont', '-ont', default='mf',
    help='GO subontology (bp, mf, cc)')
@ck.option(
    '--test-data-name', '-td', default='test', type=ck.Choice(['test', 'time']),
    help='Test data set name')
def main(data_root, ont, test_data_name):
    train_data_file = f'{data_root}/{ont}/train_data.pkl'
    valid_data_file = f'{data_root}/{ont}/valid_data.pkl'
    test_data_file = f'{data_root}/{ont}/{test_data_name}_data.pkl'
    out_file = f'{data_root}/{ont}/{test_data_name}_predictions_interpro.pkl'
    
    go_rels = Ontology(f'{data_root}/go-basic.obo', with_rels=True)
    terms_df = pd.read_pickle(f'{data_root}/{ont}/all_terms.pkl')

    terms = terms_df['gos'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}

    train_df = pd.read_pickle(train_data_file)
    valid_df = pd.read_pickle(valid_data_file)
    train_df = pd.concat([train_df, valid_df])
    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))

    test_df = pd.read_pickle(test_data_file)

    terms_set = set(terms_dict)
    # load interpro2go
    ipr2go = {}
    with open('data/interpro2go.txt') as f:
        for line in f:
            if line.startswith('!'):
                continue
            it = line.strip().split()
            ipr = it[0].split(':')[1]
            go_id = it[-1]
            if ipr not in ipr2go:
                ipr2go[ipr] = set()
            ipr2go[ipr].add(go_id)
    preds = []
    for i, row in test_df.iterrows():
        pred_scores = np.zeros(len(terms), dtype=np.float32)
        annots = set()
        for ipr in row['interpros']:
            if ipr not in ipr2go:
                continue
            for go_id in ipr2go[ipr]:
                annots |= go_rels.get_ancestors(go_id)
        for go_id in annots:
            if go_id in terms_dict:
                pred_scores[terms_dict[go_id]] = 1
        preds.append(pred_scores)
    test_df['preds'] = preds
    test_df.to_pickle(out_file)


if __name__ == '__main__':
    main()
