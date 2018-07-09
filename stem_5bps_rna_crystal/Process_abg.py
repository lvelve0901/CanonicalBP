#!/usr/bin/python

import os

inpdir = './Fragment_nts'
pdblist = []  #list to store pdb file names

#Walk data_directory and get list of pdb names
for file in os.listdir(inpdir):
    if file.endswith(".pdb"):
        pdblist.append(file)

pdblist.sort()  #sort the pdblist

for pdb in pdblist:
    ipdbf = os.path.join(inpdir,pdb)
    os.system("perl allAtomMeas_RNA.pl %s A.1:5 B.6:10 '((.((' ')).))'"%ipdbf)

