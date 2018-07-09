import os
import pandas as pd
from math import isnan
from commontool import read

df = pd.read_csv('Final_crystal.csv')
stems = df.stem[(df.cat == 'RNA') & (df.reso <= 3.0)].tolist()
stemlist = []

file = open('stemTofrag.txt','wt')
for stem in stems:
    stem = str(stem)
    if stem == 'nan':
        continue
    st = stem.split()
    for s in st:
        if int(s.split('_')[-1]) >= 3:
            stemlist.append(s)
            print(s)
            file.write('%s\n'%s)

file.close()
