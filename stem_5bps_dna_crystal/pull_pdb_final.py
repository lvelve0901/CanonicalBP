import pandas as pd

df = pd.read_table('Output_seq.txt',delim_whitespace=True)
pdblist = list(set(df.pdbid.tolist()))
pdblist.sort()

file = open('pdb_final.txt','wt')
for pdb in pdblist:
    file.write("%s\n"%pdb)

file.close()
