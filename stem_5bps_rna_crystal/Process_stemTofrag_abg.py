#!/usr/bin/python

import os
import sys
import time
from pdblib.base import *
from commontool import read, readchar


inpdir = '../Stem'
oupdir = './Fragment_abg'
ouplist = []  #list to store pdb file information
pdblist = []  #list to store pdb file names

#Walk data_directory and get list of pdb names
for file in read('stemTofrag.txt'):
    pdblist.append(file[0])

pdblist.sort()  #sort the pdblist

tic = time.time()  #timer: start

for idx, pdb in enumerate(pdblist):
    print ("--- Working on [%s] (%d of %d) ---"%(pdb,idx+1,len(pdblist)))
    pdbid = pdb.split('_')[0]
    stemid = pdb.split('_')[1] 
    ipdbf = os.path.join(inpdir,pdb+".pdb")
    mol = Mol(ipdbf)
    seg0 = mol.segs[0]
    seg1 = mol.segs[1]
    reses0 = seg0.reses
    reses1 = seg1.reses
    ch0 = seg0.chid
    ch1 = seg1.chid
    for i in range(len(reses0)-4):
        newmol = Mol()
        newmol.segs = [Segment(),Segment()]
        newseg0 = newmol.segs[0]
        newseg1 = newmol.segs[1]
        newseg0.reses = [Residue(),Residue(),Residue(),Residue(),Residue()]
        newseg1.reses = [Residue(),Residue(),Residue(),Residue(),Residue()]
        cbpnm0 = reses0[i+2].name
        cbpnm1 = reses1[-i-3].name
        cbpid0 = str(reses0[i+2].resi)
        cbpid1 = str(reses1[-i-3].resi)
        if abs(reses0[i+4].resi-reses0[i].resi) != 4 or abs(reses1[-i-1].resi-reses1[-i-5].resi) != 4:  #check point: continuous index for each strand
            print(">>>>>>>>>>>>>>>>>Abnormal index of residue!")
            print("Position: bp0 "+cbpid0+" bp1 "+cbpid1)
            #sys.exit()
        if (cbpnm0=="A" and cbpnm1=="U") or (cbpnm1=="A" and cbpnm0=="U"):
            cbpmf = "AU"
        elif (cbpnm0=="G" and cbpnm1=="C") or (cbpnm1=="G" and cbpnm0=="C"):
            cbpmf = "GC"
        elif (cbpnm0=="G" and cbpnm1=="U") or (cbpnm1=="G" and cbpnm0=="U"):
            cbpmf = "GU"
        unitid = pdbid+"_"+stemid+"_"+cbpmf+"_"+ch0+"_"+cbpnm0+"_"+cbpid0+"_"+ch1+"_"+cbpnm1+"_"+cbpid1
        restart = False
        if cbpnm0 == "A" or cbpnm0 == "G":
            for j in range(5):
                if (reses0[i+j].name == "A" and reses1[-i-j-1].name == "U") or (reses0[i+j].name == "U" and reses1[-i-j-1].name == "A") or (reses0[i+j].name == "G" and reses1[-i-j-1].name == "C") or (reses0[i+j].name == "C" and reses1[-i-j-1].name == "G") or (reses0[i+j].name == "G" and reses1[-i-j-1].name == "U") or (reses0[i+j].name == "U" and reses1[-i-j-1].name == "G"):
                    newseg0.reses[j].atoms = reses0[i+j].atoms
                    newseg0.reses[j].name = reses0[i+j].name
                    newseg0.reses[j].resi = j+1
                    newseg1.reses[j].atoms = reses1[-i+j-5].atoms
                    newseg1.reses[j].name = reses1[-i+j-5].name
                    newseg1.reses[j].resi = j+6
                else:
                    print("Unexpected base pair motif! Skip!")
                    print("Position: bp0 "+str(reses0[i+j].resi)+" "+reses0[i+j].name+" bp1 "+str(reses1[-i-j-1].resi)+" "+reses1[-i-j-1].name)
                    if reses0[i+j].name == reses1[-i-j-1].name:
                        os.system("cp %s Repair/"%ipdbf)
                    restart = True
                    break
            if restart:
                continue
            seq = [reses0[i].name[0],reses0[i+1].name[0],reses0[i+2].name[0],reses0[i+3].name[0],reses0[i+4].name[0]]
            newseg0.chid = "A"
            newseg1.chid = "B"
            newmol.write(oupdir+"/"+unitid+".pdb")
            ouplist.append([unitid,pdbid,stemid,cbpmf,''.join(seq[1:4]),''.join(seq)])
        elif cbpnm1 == "A" or cbpnm1 == "G":
            for j in range(5):
                if (reses0[i+j].name == "A" and reses1[-i-j-1].name == "U") or (reses0[i+j].name == "U" and reses1[-i-j-1].name == "A") or (reses0[i+j].name == "G" and reses1[-i-j-1].name == "C") or (reses0[i+j].name == "C" and reses1[-i-j-1].name == "G") or (reses0[i+j].name == "G" and reses1[-i-j-1].name == "U") or (reses0[i+j].name == "U" and reses1[-i-j-1].name == "G"):
                    newseg0.reses[j].atoms = reses1[-i+j-5].atoms
                    newseg0.reses[j].name = reses1[-i+j-5].name
                    newseg0.reses[j].resi = j+1
                    newseg1.reses[j].atoms = reses0[i+j].atoms
                    newseg1.reses[j].name = reses0[i+j].name
                    newseg1.reses[j].resi = j+6
                else:
                    print("Unexpected base pair motif! Skip!")
                    print("Position: bp0 "+str(reses0[i+j].resi)+" "+reses0[i+j].name+" bp1 "+str(reses1[-i-j-1].resi)+" "+reses1[-i-j-1].name)
                    if reses0[i+j].name == reses1[-i-j-1].name:
                        os.system("cp %s Repair/"%ipdbf)
                    restart = True
                    break
            if restart:
                continue
            seq = [reses1[-i-5].name[0],reses1[-i-4].name[0],reses1[-i-3].name[0],reses1[-i-2].name[0],reses1[-i-1].name[0]]
            newseg0.chid = "A"
            newseg1.chid = "B"
            newmol.write(oupdir+"/"+unitid+".pdb")
            ouplist.append([unitid,pdbid,stemid,cbpmf,''.join(seq[1:4]),''.join(seq)])

toc = time.time()

print(">>>>>> Total time used %f >>>>>>"%(toc-tic))
