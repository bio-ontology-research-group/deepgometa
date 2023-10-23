#!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from utils import Ontology, NAMESPACES

@ck.command()
@ck.option(
    '--data-root', '-d', default='data')
@ck.option('--in-file', '-i', help='Input file', required=True)
@ck.option('--out-file', '-o', help='Output file', required=True)
def main(data_root, in_file, out_file):
    go = Ontology(f'{data_root}/go.obo', with_rels=True)
    with open(in_file, 'r') as f, open(out_file, 'w') as w:
        for line in f:
            it = line.strip().split()
            go_set = set(it[1:])
            for go_id in it[1:]:
                if go_id in go_set:
                    ancestors = go.get_ancestors(go_id)
                    ancestors.discard(go_id)
                    go_set -= ancestors
                    
            w.write(it[0])
            for go_id in go_set:
                w.write(' ' + go_id)
            w.write('\n')
