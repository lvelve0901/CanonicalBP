#!/usr/bin/python

import pandas as pd

df1 = pd.read_csv('All_crystal.csv')
df2 = pd.read_csv('Stem_crystal.csv')
df3 = pd.read_csv('Bulge_crystal.csv')
df4 = pd.read_csv('Iloop_crystal.csv')
df5 = pd.read_csv('Hairpin_crystal.csv')
df6 = pd.read_csv('dna.csv')
df7 = pd.read_csv('protein#dna.csv')
df8 = pd.read_csv('rna.csv')
df9 = pd.read_csv('protein#rna.csv')
df10 = pd.read_csv('hybrid.csv')


df2 = df2.sort(['pdb']).reset_index(drop=True)
df3 = df3.sort(['pdb']).reset_index(drop=True)
df4 = df4.sort(['pdb']).reset_index(drop=True)
df5 = df5.sort(['pdb']).reset_index(drop=True)

df11 = pd.concat([df6,df7,df8,df9,df10]).sort(['pdb']).reset_index(drop=True)

df = pd.concat([df2.pdb, df1.molecule_type, df11.subcat, df1.resolution, df11.detail, df2.stem, df3.bulge, df4.iloop, df5.hairpin],axis=1)

df.to_csv('Final_crystal.csv',index=False,header = ['pdb','cat','subcat','reso','detail','stem','bulge','iloop','hairpin'])
