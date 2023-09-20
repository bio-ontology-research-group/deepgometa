#!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from collections import Counter, deque
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
    '--data-file', '-ndf', default='data/swissprot_exp_2023_03.pkl',
    help='Uniprot KB, generated with uni2pandas.py')
def main(data_root, go_file, data_file):
    go = Ontology(go_file, with_rels=True)
    logging.info('GO loaded')
    
    df = pd.read_pickle(data_file)
    
    print("DATA FILE", len(df))
    
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
            
        tdf = df.iloc[index]
        
        # Sort GO classes by their number of annotations
        terms = [g_id for g_id, c in sorted(cnt.most_common(), key=lambda x: (x[1], x[0]))]
        interpros = list(iprs.keys())

        print(f'Number of {ont} terms {len(terms)}')
        print(f'Number of {ont} iprs {len(iprs)}')

        terms_df = pd.DataFrame({'gos': terms})
        terms_df.to_pickle(out_terms_file)
        iprs_df = pd.DataFrame({'interpros': interpros})
        iprs_df.to_pickle(f'data/{ont}/interpros.pkl')
        
        # Split train/valid/test
        proteins = tdf['proteins']
        prot_set = set(proteins)
        prot_idx = {v:k for k, v in enumerate(proteins)}
        sim = {}
        train_prots = set()
        with open('data/swissprot_exp_meta.sim') as f:
            for line in f:
                it = line.strip().split('\t')
                p1, p2, score = it[0], it[1], float(it[3]) / 100.0
                if p1 == p2:
                    continue
                if p1 not in prot_set or p2 not in prot_set:
                    continue
                if p1 not in sim:
                    sim[p1] = []
                if p2 not in sim:
                    sim[p2] = []
                sim[p1].append(p2)
                sim[p2].append(p1)

        used = set()
        def dfs(prot):
            stack = deque()
            stack.append(prot)
            used.add(prot)
            prots = []
            while len(stack) > 0:
                prot = stack.pop()
                prots.append(prot)
                used.add(prot)
                if prot in sim:
                    for p in sim[prot]:
                        if p not in used:
                            used.add(p)
                            stack.append(p)
            return prots

        groups = []
        for p in proteins:
            if p not in used:
                group = dfs(p)
                groups.append(group)
        print(len(proteins), len(groups))
        index = np.arange(len(groups))
        np.random.seed(seed=0)
        np.random.shuffle(index)
        train_n = int(len(groups) * 0.9)
        valid_n = int(train_n * 0.9)
    
        train_index = []
        valid_index = []
        test_index = []
        for idx in index[:valid_n]:
            for prot in groups[idx]:
                train_index.append(prot_idx[prot])
        for idx in index[valid_n:train_n]:
            for prot in groups[idx]:
                valid_index.append(prot_idx[prot])
        for idx in index[train_n:]:
            for prot in groups[idx]:
                test_index.append(prot_idx[prot])
                
        train_index = np.array(train_index)
        valid_index = np.array(valid_index)
        test_index = np.array(test_index)

        train_df = tdf.iloc[train_index]        
        valid_df = tdf.iloc[valid_index]
        test_df = tdf.iloc[test_index]

        train_df.to_pickle(f'data/{ont}/train_data.pkl')
        valid_df.to_pickle(f'data/{ont}/valid_data.pkl')
        test_df.to_pickle(f'data/{ont}/test_data.pkl')
        

        print(f'Train/Valid/Test proteins for {ont} {len(train_df)}/{len(valid_df)}/{len(test_df)}')


if __name__ == '__main__':
    main()
