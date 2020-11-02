import pandas as pd
import sys

ml = pd.read_csv(sys.argv[1], sep='\t')
gs = set(ml.Hugo_Symbol.values)
for g in gs:
    print(g)
