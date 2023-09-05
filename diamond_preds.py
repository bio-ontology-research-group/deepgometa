#!/usr/bin/env python

import numpy as np
import pandas as pd
import click as ck
from sklearn.metrics import classification_report
from sklearn.metrics.pairwise import cosine_similarity
import sys
from collections import deque
import time
import logging
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.spatial import distance
from scipy import sparse
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
from matplotlib import pyplot as plt
import os

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@ck.command()
@ck.option(
    '--data-root', '-dr', default='data',
    help='Prediction model')
@ck.option(
    '--ont', '-ont', default='mf',
    help='Prediction model')
def main(data_root, ont):
    train_data_file = f'data/{ont}/train_data.pkl'
    valid_data_file = f'{data_root}/{ont}/valid_data.pkl'
    test_data_file = f'{data_root}/{ont}/test_data.pkl'
    terms_file = f'{data_root}/{ont}/terms.pkl'
    
    go_rels = Ontology(f'data/go-basic.obo', with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['gos'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}

    train_df = pd.read_pickle(train_data_file)
    valid_df = pd.read_pickle(valid_data_file)
    train_df = pd.concat([train_df, valid_df])
    
    test_df = pd.read_pickle(test_data_file)
    
    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    
    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i

    prot_ac_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_ac_index[row.accessions[0]] = i

    dsets = {'train': train_df, 'valid': valid_df, 'test': test_df}
    
    # BLAST Similarity (Diamond)
    for dset in ('test',):
        diamond_scores_file = f'{data_root}/{ont}/{dset}_diamond.res'
        diamond_scores = {}
        with open(diamond_scores_file) as f:
            for line in f:
                it = line.strip().split()
                if it[0] == it[1]: # Ignore same proteins (for train)
                    continue
                if it[0] not in diamond_scores:
                    diamond_scores[it[0]] = {}
                diamond_scores[it[0]][it[1]] = float(it[2])

        # scores_file = f'{data_root}/{ont}/{dset}_foldseek.m8'
        # scores = {}
        # with open(scores_file) as f:
        #     for line in f:
        #         it = line.strip().split()
        #         if it[0] == it[1]: # Ignore same proteins (for train)
        #             continue
        #         p1, p2 = it[0].split('.')[0], it[1].split('.')[0]
        #         if p1 not in scores:
        #             scores[p1] = {}
        #         scores[p1][p2] = float(it[2])

        # fold_preds = []
        # print(f'Preds for {dset}')
        # df = dsets[dset]
        # for i, row in enumerate(df.itertuples()):
        #     annots = {}
        #     prop_annots = {}
        #     prot_id = row.accessions[0]
        #     if prot_id in scores:
        #         sim_prots = scores[prot_id]
        #         allgos = set()
        #         total_score = 0.0
        #         for p_id, score in sim_prots.items():
        #             if p_id in prot_ac_index:
        #                 allgos |= annotations[prot_ac_index[p_id]]
        #                 total_score += score
        #         allgos = list(sorted(allgos))
        #         sim = np.zeros(len(allgos), dtype=np.float32)
        #         for j, go_id in enumerate(allgos):
        #             s = 0.0
        #             for p_id, score in sim_prots.items():
        #                 if p_id in prot_ac_index and go_id in annotations[prot_ac_index[p_id]]:
        #                     s += score
        #             sim[j] = s / total_score
        #         for go_id, score in zip(allgos, sim):
        #             annots[go_id] = score
        #         prop_annots = annots.copy()

        #         for go_id, score in annots.items():
        #             for sup_go in go_rels.get_ancestors(go_id):
        #                 if sup_go in prop_annots:
        #                     prop_annots[sup_go] = max(prop_annots[sup_go], score)
        #                 else:
        #                     prop_annots[sup_go] = score
        #     fold_preds.append(prop_annots)

        diam_preds = []
        print(f'Diamond preds for {dset}')
        df = dsets[dset]
        for i, row in enumerate(df.itertuples()):
            annots = {}
            prop_annots = {}
            prot_id = row.proteins
            # DiamondScore
            if prot_id in diamond_scores:
                sim_prots = diamond_scores[prot_id]
                allgos = set()
                total_score = 0.0
                for p_id, score in sim_prots.items():
                    allgos |= annotations[prot_index[p_id]]
                    total_score += score
                allgos = list(sorted(allgos))
                sim = np.zeros(len(allgos), dtype=np.float32)
                for j, go_id in enumerate(allgos):
                    s = 0.0
                    for p_id, score in sim_prots.items():
                        if go_id in annotations[prot_index[p_id]]:
                            s += score
                    sim[j] = s / total_score
                for go_id, score in zip(allgos, sim):
                    annots[go_id] = score

                prop_annots = annots.copy()
                for go_id, score in annots.items():
                    for sup_go in go_rels.get_ancestors(go_id):
                        if sup_go in prop_annots:
                            prop_annots[sup_go] = max(prop_annots[sup_go], score)
                        else:
                            prop_annots[sup_go] = score
            diam_preds.append(prop_annots)

        # new_preds = []
        # for dp, fp in zip(diam_preds, fold_preds):
        #     for go_id, score in fp.items():
        #         if go_id in dp:
        #             dp[go_id] = max(dp[go_id], score)
        #         else:
        #             dp[go_id] = score
        #     new_preds.append(dp)
        
        df['diam_preds'] = diam_preds
        # df['fold_preds'] = fold_preds
        # df['diam_fold_preds'] = new_preds
        filename = f'{data_root}/{ont}/{dset}_data_diam.pkl'
        df.to_pickle(filename)

if __name__ == '__main__':
    main()
