"""TumOnc package

Fedor.Bezrukov@gmail.com
"""

import pandas as pd
import numpy as np

from TumOnc.critstat import critstat
from TumOnc.context import std_ctx3, std_ctx5, snp_id
from TumOnc.genome import read_genes_default
from TumOnc.binning import binned_ns, binned_ns_load

__version__ = "0.01"

# __all__ = ["critstat","context","genome","pmodel","binning"]


# All the work

# Read and standartize mutations


def read_mutations(file="BCC_MUT_combined_context.maf"):
    muts = pd.read_csv(file, sep='\t', low_memory=False)
    muts.rename(columns={'Hugo_Symbol': 'gene',
                         'Tumor_Sample_Barcode': 'patient'},
                inplace=True)
    muts['ctx3'] = muts['ref_context'].str.upper().apply(lambda x: std_ctx3[x[10-1:10+2]])
    muts['ctx5'] = muts['ref_context'].str.upper().apply(lambda x: std_ctx5[x[10-2:10+3]])
    # Not using categories for now
    # for c in ['gene','patient','Variant_Type','Variant_Classification','ctx3','ctx5']:
    #     muts[c] = muts[c].astype("category")
    # Cleanup
    for c in range(1, 11):
        muts.drop(['nt'+str(c)+'m', 'nt'+str(c)+'p'], axis=1)
    muts.drop('nt00', axis=1)
    muts['snp_id']=pd.Series(snp_id)[
        muts.apply(lambda x:x['Reference_Allele']+str(x['Tumor_Seq_Allele2']),axis=1)
    ].values
    return muts


def get_snp_mutations(bcc):
    return bcc[bcc.Variant_Type == 'SNP']


def get_missense_mutations(bcc):
    return bcc[bcc.Variant_Classification == 'Missense_Mutation']


def get_silent_mutations(bcc):
    return bcc[bcc.Variant_Classification == 'Silent']


def clean_mutations(bcc, genes):
    # Filter out "good" genes
    bcc = bcc[bcc.apply(lambda x:x['gene'] in genes, axis=1)].copy()
    # Recategorize to drop empty genes and patients
    # bcc['gene'] = bcc['gene'].astype('str').astype('category')
    # bcc['patient'] = bcc['patient'].astype('str').astype('category')
    # Drop Unknown mutations
    bcc.dropna(subset=['Tumor_Seq_Allele2'], inplace=True)
    return bcc


def search_oncogenes(nobin, pmodel_nobin=False,
                     pval4=False, minmuts=3):
    if not pmodel_nobin:
        pmodel_nobin = nobin
    bcc = nobin.muts
    bcc['ppos'] = bcc.apply(pmodel_nobin.ppos, axis=1)

    def mhit_pvalue(g, patient_factor, patient_id):
        plist = patient_id[g['patient']].values
        pval = g.iloc[0].at['ppos']
        cs = critstat(pval, patient_factor.values)
        pguess = cs.guessp(plist)
        if pval4:
            return cs.pvalue4(plist), pguess
        else:
            return cs.pvalue3(plist), pguess
        
    hits = bcc.groupby( ('Chromosome','Start_position') )
    pvals = {}
    guesspvals = {}
    lastchr = ''
    for h, g in hits:
        if len(g) >= minmuts:
            if h[0] != lastchr:
                print(h[0], end=' ', flush=True)
                lastchr = h[0]
            # print(h,end=' -> ')
            pvals[h], guesspvals[h] = mhit_pvalue(g, nobin.patient_factor,
                                                  nobin.patient_id)
            # print(pvals[h])
        
    pvals_s = pd.DataFrame({'pval': pd.Series(pvals),
                            'guesspval': pd.Series(guesspvals)})
    pvals_s['nhits'] = hits.size()
    pvals_s = pvals_s.sort('pval')
    pvals_s['qval']=pvals_s.pval.values*np.sum(nobin.N)/range(1,len(pvals_s)+1)
    pvals_s['mut_id'] = None
    for k in pvals:
        pvals_s.loc[k, 'mut_id'] = hits.groups[k][0]
    pvals_s['ctx5'] = bcc.loc[pvals_s.mut_id, 'ctx5'].values
    pvals_s['gene'] = bcc.loc[pvals_s.mut_id, 'gene'].values
    pvals_s['cDNA_Change']=bcc.loc[pvals_s.mut_id,'cDNA_Change'].values
    pvals_s['ppos'] = bcc.loc[pvals_s.mut_id, 'ppos'].values
    if hasattr(pmodel_nobin.pmodel,'expcv'):
        pvals_s['cov_fact'] = pvals_s['gene'].map(nobin.bin_gene).map(pd.Series(pmodel_nobin.pmodel.expcv))
    return pvals_s,hits

