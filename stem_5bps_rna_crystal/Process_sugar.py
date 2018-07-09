#!/usr/bin/python

import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import re
import learnna_json as lna_json
from pdblib.num import *

distlist_Pu = [["C1'","N9"],["C3'","O3'"],["C4'","C5'"],["C1'","C2'"],["C2'","C3'"],["C3'","C4'"],["C4'","O4'"]]
distlist_Py = [["C1'","N1"],["C3'","O3'"],["C4'","C5'"],["C1'","C2'"],["C2'","C3'"],["C3'","C4'"],["C4'","O4'"]]
distlist_Bp = [["C1'","C1'"]]
anglist_Pu = [["C4'","C1'","N9"],["O4'","C1'","N9"],["C2'","C1'","N9"],["C2'","C3'","O3'"],["C4'","C3'","O3'"],["C3'","C4'","C5'"],["O4'","C4'","C5'"]]
anglist_Py = [["C4'","C1'","N1"],["O4'","C1'","N1"],["C2'","C1'","N1"],["C2'","C3'","O3'"],["C4'","C3'","O3'"],["C3'","C4'","C5'"],["O4'","C4'","C5'"]]

pdb_dir = './Fragment_nts'
oupname = 'Output_sugar'
pdblist = []  #list to store PDB file names
output = []  #list to store output information

#Walk data_directory and get list of pdb names
for file in os.listdir(pdb_dir):
	if file.endswith(".pdb"):
		pdblist.append(file)

pdblist.sort()

#Loop over list of files in data dir, process with DSSR to get json file
for idx, pdb in enumerate(pdblist):
    print(">>>>>>Working on %s [%d/%d]>>>>>>"%(pdb,idx,len(pdblist)))
    unitid = str(pdb)[:-4]
    confAB = []
    confAB.append(unitid)
    mol = Mol('%s/%s'%(pdb_dir,pdb))
    resA = mol.segs[0].reses[2]
    resB = mol.segs[1].reses[2]
    pairAB = [resA, resB]
    for dd in distlist_Pu:
        at0 = getat(mol,resA.resi,dd[0])
        at1 = getat(mol,resA.resi,dd[1])
        if at0 is not None and at1 is not None:
            confAB.append(dist(at0,at1))
        else:
            print("===================")
            print("No such atom: %s %s"%(dd[0],dd[1]))
            print("===================")
            confAB.append(None)
    for aa in anglist_Pu:
        at0 = getat(mol,resA.resi,aa[0])
        at1 = getat(mol,resA.resi,aa[1])
        at2 = getat(mol,resA.resi,aa[2])
        if at0 is not None and at1 is not None and at2 is not None:
           confAB.append(angle(at0,at1,at2))
        else:
           print("==========================")
           print(">>> No such atom: %s %s %s"%(aa[0],aa[1],aa[2]))
           print("==========================")
           confAB.append(None)
    for dd in distlist_Py:
        at0 = getat(mol,resB.resi,dd[0])
        at1 = getat(mol,resB.resi,dd[1])
        if at0 is not None and at1 is not None:
            confAB.append(dist(at0,at1))
        else:
            print("===================")
            print("No such atom: %s %s"%(dd[0],dd[1]))
            print("===================")
            confAB.append(None)
    for aa in anglist_Py:
        at0 = getat(mol,resB.resi,aa[0])
        at1 = getat(mol,resB.resi,aa[1])
        at2 = getat(mol,resB.resi,aa[2])
        if at0 is not None and at1 is not None and at2 is not None:
            confAB.append(angle(at0,at1,at2))
        else:
            print("==========================")
            print(">>> No such atom: %s %s %s"%(aa[0],aa[1],aa[2]))
            print("==========================")
            confAB.append(None)
    for dd in distlist_Bp:
        at0 = getat(mol,resA.resi,dd[0])
        at1 = getat(mol,resB.resi,dd[1])
        if at0 is not None and at1 is not None:
            confAB.append(dist(at0,at1))
        else:
            print("===================")
            print("No such atom: %s %s"%(dd[0],dd[1]))
            print("===================")
            confAB.append(None)
    output.append(confAB)
    

file_nts = open('%s.txt'%oupname,'w')

file_nts.write("unitid pur_C1N9 pur_C3O3 pur_C4C5 pur_C1C2 pur_C2C3 pur_C3C4 pur_C4O4 pur_C4C1N9 pur_O4C1N9 pur_C2C1N9 pur_C2C3O3 pur_C4C3O3 pur_C3C4C5 pur_O4C4C5 pym_C1N1 pym_C3O3 pym_C4C5 pym_C1C2 pym_C2C3 pym_C3C4 pym_C4O4 pym_C4C1N1 pym_O4C1N1 pym_C2C1N1 pym_C2C3O3 pym_C4C3O3 pym_C3C4C5 pym_O4C4C5 bp_c1pc1p\n")
for nt in output:
	file_nts.write('%30s %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %5.3f\n'%(nt[0],nt[1],nt[2],nt[3],nt[4],nt[5],nt[6],nt[7],nt[8],nt[9],nt[10],nt[11],nt[12],nt[13],nt[14],nt[15],nt[16],nt[17],nt[18],nt[19],nt[20],nt[21],nt[22],nt[23],nt[24],nt[25],nt[26],nt[27],nt[28],nt[29]))

file_nts.close()
