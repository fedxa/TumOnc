from __future__ import print_function
from __future__ import division

import numpy as np
import math as ma

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation

import csv
# import gzip
import os
import pickle

from TumOnc.context import std_ctx3, std_ctx5

########################################################################
# Genome initialization


def read_genome(file="human_g1k_v37.fasta", full=False):
    """Read the genome sequence

    Returns the genome as a dictionary of SeqRecord objects labeled by
    the chromosome name -- in particular, hg['1'].id == '1'

    if full=True is specified the file is really read, otherwise it is read
    on demand for each call
    """
    if full:
        # Internal gzip implementaton fails on python3?!
        # with gzip.open(file,"rb") as fd:
        if file.endswith("gz"):
            with os.popen("zcat "+file, 'r') as fd:
                hg = SeqIO.to_dict(SeqIO.parse(fd, "fasta", generic_dna))
        else:
            with open(file, 'r') as fd:
                hg = SeqIO.to_dict(SeqIO.parse(fd, "fasta", generic_dna))
    else:
        hg = SeqIO.index_db(file+".idx", file, "fasta", alphabet=generic_dna)
    return hg


class gene:

    # def __str__(self):
    #     return self.__str__()
    def __repr__(self):
        return """Gene
   chr: '%s'
   loc: '%s'
   covariates :
         lexpr ='%s'
         rept ='%s'
         hic  ='%s'
        """ % (getattr(self, 'chr', '?'),
               getattr(self, 'loc', '?'),
               getattr(self, 'lexpr', 'nan'),
               getattr(self, 'rept', 'nan'),
               getattr(self, 'hic', 'nan'))

    def __init__(self):
        self.chr = ''
        self.lexpr = np.nan
        self.rept = np.nan
        self.hic = np.nan

    def fill_covariates(self, e, r, h):
        try:
            self.lexpr = ma.log(float(e))
        except:
            self.lexpr = np.nan
        try:
            self.rept = float(r)
        except:
            self.rept = np.nan
        try:      
            self.hic = float(h)
        except:
            self.hic = np.nan
            
    def good(self):
        return (self.chr and
                hasattr(self, 'loc') and
                not np.isnan(self.lexpr) and
                not np.isnan(self.rept) and
                not np.isnan(self.hic))

    def extract(self, hg):
        return self.loc.extract(hg[self.chr])

    # Returned scaled values of the covariates
    def Vlexpr(self, scale=20):
        return self.lexpr

    def Vrept(self, scale=2000):
        return self.rept

    def Vhic(self, scale=100):
        return self.hic


# What is the meaning of interleaved genes?
def read_gene_exones(file="Agil_RS_names.bed"):
    """Reads genes as features with locations for the whole genome"""
    with open(file, "r") as fd:
        lg = sorted(csv.reader(fd, dialect='excel-tab'),
                    key=lambda x: '{} {} {:010d}'.format(x[0], x[3],
                                                         int(x[1])))
        # Sort as chr:name:position
    genes = {}
    for chr, spfield, epfield, name in lg:
        sp = int(spfield)+2  # Ok, these are based on some strange asumptions
        ep = int(epfield)-2  # convert _exactly_ to exon locations (?)
        if name in genes:
            g = genes[name]
        else:
            g = gene()  # Create new gene entry
            g.chr = chr
            genes[name] = g
        if hasattr(g, 'loc'):
            g.loc += FeatureLocation(sp, ep)
        else:
            g.loc = FeatureLocation(sp, ep)
    return genes

def fill_gene_covariates(genes,file="gene.covariates.txt"):
    with open(file,"r") as fd:
        for g in csv.DictReader(fd,dialect='excel-tab'):
            gn = g['gene']
            if gn not in genes:
                genes[gn] = gene()
            genes[gn].fill_covariates(g['expr'],g['reptime'],g['hic'])


def normalize_gene_covariates(genes):
    V = np.empty((3, len(genes)))
    i = 0
    for gn in genes:
        V[0, i] = genes[gn].lexpr
        V[1, i] = genes[gn].rept
        V[2, i] = genes[gn].hic
        i += 1
    for i in range(3):
        V[i,:] -= np.nanmean(V[i,:])
        V[i,:] /= np.nanstd(V[i,:], ddof=1)
    i = 0
    for gn in genes:
        genes[gn].lexpr = V[0, i]
        genes[gn].rept = V[1, i]
        genes[gn].hic = V[2, i]
        i += 1
    return genes


def count_contexts_in_exon(counts,loc,chrseq,ctxstd,ctxlen):
    ct2 = ctxlen//2
    start = loc.start-ct2
    end = loc.end+ct2
    assert start>=0 and end<=len(chrseq)
    exonstr = str(chrseq[start:end].seq)
    for pos in range(ct2,len(exonstr)-ct2):
        counts[ctxstd[exonstr[pos-ct2:pos+ct2+1]]] += 1

def count_contexts_in_gene_cds(loc,chrseq,ctxstd,ctxlen):
    counts = { k : 0 for k in ctxstd.values() }
    for p in loc.parts:
        count_contexts_in_exon(counts,p,chrseq,ctxstd,ctxlen)
    return counts

def fill_gene_contexts(genes,ctxstd,ctxlen,ctxattrname,
                       genome,
                       prefix=''):
    fn = prefix+'.'+ctxattrname+'.pcl'
    try: # Read existing data if present
        with open(fn,'rb') as f:
            ctxcounts = pickle.load(f)
        print("Read precalculated counts from "+fn)
        for g in ctxcounts:
            setattr(genes[g],ctxattrname,ctxcounts[g])
    except: # Otherwise recalculate
        hg = read_genome(genome)
        ctxcounts = {}
        for chrid in hg.keys():
            gl = [g for g in genes if genes[g].chr==chrid]
            if gl: # If there are genes in the current chromosome
                print("Counting in chr",chrid)
                chrseq=hg[chrid]
                for g in gl:
                    # print(" ",g)
                    if hasattr(genes[g],'loc'): # For genes with locations
                        counts = count_contexts_in_gene_cds(genes[g].loc,chrseq,
                                                            ctxstd,ctxlen)
                        setattr(genes[g],ctxattrname,counts)
                        ctxcounts[g] = counts
        with open(fn, 'wb') as f:  # And save for future
            print('Saving context counts to '+fn)
            pickle.dump(ctxcounts, f)


def read_genes_default(exones="Agil_RS_names.bed",
                       covariates="gene.covariates.txt",
                       genome="human_g1k_v37.fasta",
                       normalize=True):
    """Read all the standard data about genes
    _and_ cleans them by dropping genes with unknown covariates!
    """
    genes = read_gene_exones(exones)
    fill_gene_covariates(genes, covariates)
    if normalize:
        normalize_gene_covariates(genes)
    fill_gene_contexts(genes, std_ctx3, 3, 'ctx3counts', genome, prefix=exones)
    fill_gene_contexts(genes, std_ctx5, 5, 'ctx5counts', genome, prefix=exones)
    for g in genes:
        if hasattr(genes[g], 'loc'):
            genes[g].len = len(genes[g].loc)
    for g in list(genes.keys()):
        if not genes[g].good() or genes[g].lexpr == 0 or genes[g].rept == 0:
            del genes[g]
    return genes


def genes_def_cat(genes):
    """Returns the 'standard' numbering for the genes"""
    return {g: i for i, g in enumerate(genes.keys())}
