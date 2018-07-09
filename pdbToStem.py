#!/usr/bin/python

import os
import sys
import time
import pandas as pd
from pdblib.base import *
from commontool import read, readchar

def addTER(filename):  #Add TER in the dssr-stem file 
    strings = read(filename)  #inport stem file into string list
    chars = readchar(filename)  #import stem file into char list
    ters = []  #line index that should add ter
    models = []  #number of MODEL line
    bpnums = [] #number of bps in each MODEL
    count = 0  #initialize count of residue id
    bpnum = -1  #initialize num of bps
    for i in range(len(strings)):
        if strings[i][0] == 'REMARK':
            bpnum = int(strings[i][2][4:])  #read num of bps
            bpnums.append(bpnum)
            models.append(strings[i-1][1])  #read model index
            count = 0  #reinitialize count of residue id
        if strings[i][0] in ['ATOM','HETATM'] and strings[i+1][0] in ['ATOM','HETATM']:
            if chars[i][17:27] != chars[i+1][17:27]:  #read column of resi
                count += 1  #count + 1 if the resi change
                if count == bpnum:  #if count of resi = num of bp
                    ters.append(i)
    for i in range(len(ters)):
        os.system("sed -i '%d i\TER' dssr-stems.pdb"%(ters[i]+i+2))
    return ters, bpnums, models

inpdir = './Crystal'
jsondir = './Json'
oupdir = './Stem'
curdir = os.getcwd()
stemdic = {}  #dictionary to store pdb -- stem information
pdblist = []  #list to store pdb file names
jsonlist = []  #list to store existing json file names

#Walk pdb_directory and get list of pdb names
for file in read('Pdbinfo/All_crystal.txt'):
    pdblist.append(file[0])

pdblist.sort()  #sort the pdblist

#Walk json_directory to check if the pdb exist in the json file
for file in os.listdir(jsondir):
    jsonlist.append(file)

jsonlist.sort()  #sort the pdblist

tic = time.time()  #timer: start

for idx, pdbid in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdbid,idx+1,len(pdblist)))
    pdb = pdbid + ".pdb"
    stemdic[pdbid] = ''
    #Step1: Create output json file and dssr-stems file
    ipdbf = os.path.join(inpdir,pdb)
    jsonf = pdb.replace(".pdb",".json")
    if jsonf in jsonlist:  #check if we already process the pdbfile
        continue
    os.system("x3dna-dssr --more --json --symm -i=%s -o=%s"%(ipdbf,jsonf))
    #Step2: Check if there is output dssr-stem file
    if "dssr-stems.pdb" not in os.listdir(curdir):
        os.system("cp %s %s"%(jsonf,os.path.join(jsondir,jsonf)))
        os.system('x3dna-dssr --cleanup')
        os.system("rm %s"%jsonf)
        continue
    #Step3: Add TER between two strands for each stem in stem file
    ter, bpnums, model = addTER("dssr-stems.pdb")
    if len(ter) != len(model) or len(ter) != len(bpnums):  #check point: ter bpnums and model has same number
        print(">>>Abort! Inconsistent model, ter and bps number: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
        sys.exit()
    #Step4: Split stem file into each stem using pdblib
    stems = Pdb("dssr-stems.pdb")
    for i, md in enumerate(stems.mds):
        if len(md.segs) > 2:
            print(">>>Abort! There are more than 2 chains: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        if len(md.segs[0].reses) != len(md.segs[1].reses) or len(md.segs[0].reses) != int(bpnums[i]):
            print(">>>In model %d, seg0 length: %d seg1 length: %d bpnums: %d"%(i,len(md.segs[0].reses),len(md.segs[1].reses),int(bpnums[i])))
            print(">>>Abort! Wrong position of TER: %s [%d/%d]"%(pdb,idx+1,len(pdblist)))
            sys.exit()
        opdbf = os.path.join(oupdir,pdbid+"_"+str(i)+"_"+str(len(md.segs[0].reses))+".pdb")  #pdbid-stemid-stemlen.pdb
        stemdic[pdbid] += " "+pdbid+"_"+str(i)+"_"+str(len(md.segs[0].reses))
        md.write("%s"%opdbf)
    os.system("cp %s %s"%(jsonf,os.path.join(jsondir,jsonf)))
    os.system("x3dna-dssr --cleanup")
    os.system("rm %s"%jsonf)
    print(stemdic[pdbid])

toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))

stem_df = pd.DataFrame()
stem_df['pdb'] = stemdic.keys()
stem_df['stem'] = stemdic.values()
stem_df.to_csv("Stem_crystal.csv",index=False)
