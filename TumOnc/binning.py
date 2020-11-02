from __future__ import print_function
from __future__ import division

import numpy as np
import pandas as pd
import pickle

from TumOnc.pmodel import pmodel3
from TumOnc.genome import genes_def_cat
from TumOnc.context import ctx5_ctx3id


# Binning

# Fill category based arrays

def fill_N(genes, ctx_cat, gene_cat):
    N = np.zeros((len(set(ctx_cat.values())), len(set(gene_cat.values()))))
    for g in genes:
        for ctx, num in genes[g].ctx5counts.items():
            N[ctx_cat[ctx], gene_cat[g]] += num
    return N

 
def fill_ns(muts,ctx_cat,gene_cat,patient_id):
    ns = np.zeros( ( 3, len(set(ctx_cat.values())), len(set(gene_cat.values())),
                     patient_id.size ) )
    for mid,mut in muts.iterrows():
        ns[ mut['snp_id'], ctx_cat[mut['ctx5']],
            gene_cat[mut['gene']], patient_id[mut['patient']] ] += 1
    return ns


def fill_Vlexpr(genes,gene_cat,scale=20):
    l = len(set(gene_cat.values()))
    vals = np.zeros(l)
    counts = np.zeros(l)
    for g in genes:
        vals[gene_cat[g]] += genes[g].Vlexpr(scale)
        counts[gene_cat[g]] += 1
    return vals/counts


def fill_Vrept(genes,gene_cat,scale=2000):
    l = len(set(gene_cat.values()))
    vals = np.zeros(l)
    counts = np.zeros(l)
    for g in genes:
        vals[gene_cat[g]] += genes[g].Vrept(scale)
        counts[gene_cat[g]] += 1
    return vals/counts


def fill_Vhic(genes,gene_cat,scale=100):
    l = len(set(gene_cat.values()))
    vals = np.zeros(l)
    counts = np.zeros(l)
    for g in genes:
        vals[gene_cat[g]] += genes[g].Vhic(scale)
        counts[gene_cat[g]] += 1
    return vals/counts


def binned_ns_load(file):
    with open(file, 'rb') as f:
        bn = pickle.load(f)
    return bn


class binned_ns:
    def __init__(self, muts, genes):
        self.muts = muts
        self.genes = genes
        self.C = len(set(muts.patient.values))
        self.patient_factor = muts.groupby('patient').size()
        self.patient_factor.sort_values(ascending=False)
        self.patient_factor = self.patient_factor/self.patient_factor.mean()
        self.patient_id = pd.Series(range(self.patient_factor.size),
                                    index=self.patient_factor.index)

    def save(self, file):
        with open(file, 'wb') as f:  # And save for future
            pickle.dump(self, f)

    def bin(self, bin_cat=None, bin_gene=None, pmodel=None):
        if bin_cat:
            self.bin_cat = bin_cat
        else:
            self.bin_cat = ctx5_ctx3id
        if bin_gene:
            self.bin_gene = bin_gene
        else:
            self.bin_gene = genes_def_cat(self.genes)
        self.ns = fill_ns(self.muts, self.bin_cat, self.bin_gene, self.patient_id)
        self.N = fill_N(self.genes, self.bin_cat, self.bin_gene)
        self.Vlexpr = fill_Vlexpr(self.genes, self.bin_gene)
        self.Vrept = fill_Vrept(self.genes, self.bin_gene)
        self.Vhic = fill_Vhic(self.genes, self.bin_gene)
        if pmodel:
            self.pmodel = pmodel
        else:
            self.pmodel = pmodel3()
        if pmodel != 0:
            self.pmodel.fit(self.ns, self.N, [self.Vlexpr, self.Vrept, self.Vhic])
        return self

    def ppos(self,mut):
        return self.pmodel.pmatrix[ self.bin_cat[mut.at['ctx5']],
                                    self.bin_gene[mut.at['gene']] ]
