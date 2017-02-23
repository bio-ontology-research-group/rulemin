#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import

import click as ck
import pandas as pd
import numpy as np
import gzip
import os
from utils import (
    EXP_CODES, get_ontology, get_anchestors, get_specifics, get_subset,
    MOLECULAR_FUNCTION, BIOLOGICAL_PROCESS, CELLULAR_COMPONENT)
from collections import deque, Counter
from multiprocessing import Pool


@ck.command()
def main():
    global go
    go = get_ontology('data/go.obo')
    global hp
    hp = get_ontology('data/hp.obo')
    run()


def run():
    functions, phenotypes = load_data()
    terms = list()
    n = len(functions)
    global counter
    counter = Counter()
    global tree
    tree = dict()
    e = 100
    term_index = dict()
    term_list = list()
    for go_id in go:
        term_index[go_id] = len(term_index)
        term_list.append(go_id)
    for hp_id in hp:
        term_index[hp_id] = len(term_index)
        term_list.append(hp_id)
    for i in xrange(n):
        funcs = set(map(lambda x: term_index[x], functions[i]))
        phenos = set(map(lambda x: term_index[x], phenotypes[i]))
        terms.append(funcs | phenos)
        for func in funcs:
            for pheno in phenos:
                counter[frozenset([func, pheno])] += 1
    for s, c in counter.items():
        if c < e:
            del counter[s]
    for s, c in counter.items():
        for term in s:
            if term_list[term] in go:
                tree[term] = set(map(
                    lambda x: term_index[x],
                    get_anchestors(go, term_list[term])))
                tree[term] |= set(map(
                    lambda x: term_index[x],
                    get_subset(go, term_list[term])))
            else:
                tree[term] = set(map(
                    lambda x: term_index[x],
                    get_anchestors(hp, term_list[term])))
                tree[term] |= set(map(
                    lambda x: term_index[x],
                    get_subset(hp, term_list[term])))

    print(len(counter))
    pool = Pool(48)
    gf = gzip.open('data/results.gz', 'w')
    while len(counter) > 0:
        cnts = pool.map(next_level, terms)
        cnt = sum(cnts)
        print(counter.most_common(10))
        print(cnt.most_common(10))
        for s, c in cnt.items():
            if c < e:
                del cnt[s]
            else:
                gf.write(c)
                for term in s:
                    gf.write('\t' + term_list[term])
                gf.write('\n')
        counter = cnt


def next_level(terms):
    cnt = Counter()
    for s in counter:
        if s.issubset(terms):
            term_set = set()
            for term in s:
                term_set |= tree[term]
            for term in terms:
                if term not in term_set:
                    ss = set(s)
                    ss.add(term)
                    cnt[frozenset(ss)] += 1
    return cnt


def uni2gene():
    df = pd.read_pickle('data/idmapping.9606.pkl')
    mapping = dict()
    for i, row in df.iterrows():
        if isinstance(row['genes'], str):
            mapping[row['accessions']] = row['genes']
    return mapping


def load_data():
    df = get_data()
    functions = list()
    phenotypes = list()
    n = len(df)
    train_n = int(0.8 * n)
    index = np.arange(n)
    np.random.seed(seed=0)
    np.random.shuffle(index)
    train_df = df.loc[index[:train_n]]
    test_df = df.loc[index[train_n:]]
    test_df.to_pickle('data/test_data.pkl')
    for i, row in train_df.iterrows():
        funcs = set()
        phenos = set()
        for func in row['functions']:
            funcs |= get_anchestors(go, func)
        for pheno in row['phenotypes']:
            phenos |= get_anchestors(hp, pheno)
        phenos.discard('HP:0000001')
        funcs.discard(MOLECULAR_FUNCTION)
        funcs.discard(BIOLOGICAL_PROCESS)
        funcs.discard(CELLULAR_COMPONENT)
        functions.append(funcs)
        phenotypes.append(phenos)
    return functions, phenotypes


def get_data():
    if os.path.exists('data/data.pkl'):
        return pd.read_pickle('data/data.pkl')
    func_annots = dict()
    pheno_annots = dict()
    mapping = uni2gene()
    with gzip.open('data/goa_human.gaf.gz') as f:
        for line in f:
            if line.startswith('!'):
                continue
            items = line.strip().split('\t')
            if items[6] not in EXP_CODES:
                continue
            prot_id = items[1]
            if prot_id not in mapping:
                print('Missing protein')
                continue
            gene_id = mapping[prot_id]
            go_id = items[4]
            if go_id not in go:
                continue
            if gene_id not in func_annots:
                func_annots[gene_id] = set()
            func_annots[gene_id].add(go_id)

    with open('data/genes_to_phenotype.txt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.strip().split('\t')
            gene_id = items[0]
            hp_id = items[3]
            if hp_id not in hp:
                continue
            if gene_id not in pheno_annots:
                pheno_annots[gene_id] = set()
            pheno_annots[gene_id].add(hp_id)

    genes = list(set(pheno_annots).intersection(set(func_annots)))
    functions = list()
    phenotypes = list()
    for gene_id in genes:
        functions.append(get_specifics(go, func_annots[gene_id]))
        phenotypes.append(get_specifics(hp, pheno_annots[gene_id]))
    df = pd.DataFrame(
        {'genes': genes, 'functions': functions, 'phenotypes': phenotypes})
    df.to_pickle('data/data.pkl')
    return df


if __name__ == '__main__':
    main()
