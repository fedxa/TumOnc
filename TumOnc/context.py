"""Various definitions related to the context handling


"""

from Bio.Seq import Seq
import pandas as pd

########################################################################
# Standard contexts and categories


def _std_context_fun(s):
    """Returns the context mapped always to the CT containing strand"""
    if s[len(s)//2] in 'CT':
        return s
    else:
        return str(Seq(s).reverse_complement())

# Standard context has 'C' or 'T' in the middle.
# If not, the reverse_complement should be used

nucl_id = {'C':0,'T':1,'A':2,'G':3}
snp_id = { 'CT':0,
           'TC':0,
           'CG':1,
           'CA':2,
           'TA':1,
           'TG':2 }
for snp in list(snp_id.keys()): # Fill up complementary cases
    snp_id[str(Seq(snp).complement())] = snp_id[snp]

# Dictionary which returns the standard length3 context
std_ctx3 = { a+b+c : _std_context_fun(a+b+c)
             for a in 'GCAT' for b in 'GCAT' for c in 'GCAT' }
def _ctx3_id_fun(c):
    return (nucl_id[c[1]]*4+nucl_id[c[0]])*4+nucl_id[c[2]]
ctx3_id = { c:_ctx3_id_fun(c) for c in set(std_ctx3.values()) }
id_ctx3 = { _ctx3_id_fun(c):c for c in set(std_ctx3.values()) }

# Dictionary which returns the standard length5 context
std_ctx5 = { a+b+c+d+e : _std_context_fun(a+b+c+d+e)
             for a in 'GCAT' for b in 'GCAT'
             for c in 'GCAT' for d in 'GCAT' for e in 'GCAT' }
def _ctx5_id_fun(c):
    return ( (((nucl_id[c[2]]*4+nucl_id[c[1]])*4+nucl_id[c[3]])*4
              +nucl_id[c[0]])*4+nucl_id[c[4]] )
ctx5_id = { c:_ctx5_id_fun(c) for c in set(std_ctx5.values()) }
id_ctx5 = { _ctx5_id_fun(c):c for c in set(std_ctx5.values()) }

ctx5_ctx3id = {ctx:ctx3_id[ctx[1:4]] for ctx in ctx5_id}



# Recheck reindexing!
def get_ctx_bins(bcc):
    """
    Doc
    """
    c3 = pd.DataFrame(bcc.groupby('ctx3').size(),columns=['ctx3count'])
    c3 = c3.reindex(ctx3_id.keys(),fill_value=0)
    c3.sort_values(by='ctx3count',ascending=False,inplace=True)
    c3['c3id'] = range(len(c3))
    c3['c3id_15'] = c3['c3id']
    c3['c3id_5'] = c3['c3id']
    c3.loc[5:,'c3id_5'] = 5
    c3.loc[15:,'c3id_15'] = 15

    c5 = pd.DataFrame(bcc.groupby('ctx5').size(),columns=['ctx5count'])
    c5 = c5.reindex(ctx5_id.keys(),fill_value=0)
    c5.sort_values(by='ctx5count',ascending=False,inplace=True)
    c5['c5id'] = range(len(c5))
    c5['c5id_15'] = c5['c5id']
    c5['c5id_5'] = c5['c5id']
    c5.loc[5:,'c5id_5'] = 5
    c5.loc[15:,'c5id_15'] = 15
    c5['ctx5'] = c5.index
    c5['ctx3'] = c5['ctx5'].str[1:4]
    c5 = c5.join(c3,on='ctx3')
    del(c5['ctx3'])
    del(c5['ctx5'])
    return c5
