#!/usr/bin/python

import json
import os
import sys
import time
from pdblib.base import *
import learnna_json as lna_json
from commontool import read, readchar


jsondir = '../Json'
stemdir = '../Stem'
oupdir = './Fragment_nts'
ouplist = []  #list to store output 5bps file information
stemlist = []  #list to store stem pdb file names

#Walk data_directory and get list of pdb names
for file in read('stemTofrag.txt'):
    stemlist.append(file[0])

stemlist.sort()  #sort the stemlist

tic = time.time()  #timer: start

for idx, stem in enumerate(stemlist):
    print ("--- Working on [%s] (%d of %d) ---"%(stem,idx+1,len(stemlist)))
    pdbid = stem.split('_')[0]
    stemid = stem.split('_')[1] 
    stemid_int = int(stemid)
    ijsonf = os.path.join(jsondir,pdbid+".json")
    istemf = os.path.join(stemdir,stem+".pdb")
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(ijsonf) as json_data: data = json.load(json_data) #read each json
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    mol = Mol(istemf)
    seg0 = mol.segs[0]
    seg1 = mol.segs[1]
    reses0 = seg0.reses
    reses1 = seg1.reses
    ch0 = seg0.chid
    ch1 = seg1.chid
    sform = "B"
    hform = najson.json_file['stems'][stemid_int]['helix_form']
    if "A" in hform: sform = "A"
    elif "Z" in hform: sform = "Z"
    elif "A" not in hform and "B" not in hform and "Z" not in hform: sform = "X"
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
        #check point: continuous index for each strand
        if abs(reses0[i+4].resi-reses0[i].resi) != 4 or abs(reses1[-i-1].resi-reses1[-i-5].resi) != 4:
            print("==========================")
            print("Abnormal index of residue!")
            print("==========================")
            print("Position: bp0 "+cbpid0+" bp1 "+cbpid1)
        #check point: center bp is canonical bp
        #make unitid
        if (cbpnm0=="DA" and cbpnm1=="DT") or (cbpnm1=="DA" and cbpnm0=="DT"): cbpmf = "AT"
        elif (cbpnm0=="DG" and cbpnm1=="DC") or (cbpnm1=="DG" and cbpnm0=="DC"): cbpmf = "GC"
        else:
            print("========================================")
            print("Unexpected center base pair motif! Skip!")
            print("========================================")
            print("Position: bp0 "+cbpid0+" "+cbpnm0+" bp1 "+cbpid1+" "+cbpnm1)
            continue
        unitid = pdbid+"_"+stemid+"_"+cbpmf+"_"+ch0+"_"+cbpnm0+"_"+cbpid0+"_"+ch1+"_"+cbpnm1+"_"+cbpid1
        #check point: swap purine and pyrimidine
        #make seq: center purine in the 1st strand
        if cbpnm0 == "DA" or cbpnm0 == "DG":
            swap = False
            seq = [reses0[i].name[1],reses0[i+1].name[1],reses0[i+2].name[1],reses0[i+3].name[1],reses0[i+4].name[1]]
        elif cbpnm1 == "DA" or cbpnm1 == "DG":
            swap = True
            seq = [reses1[-i-5].name[1],reses1[-i-4].name[1],reses1[-i-3].name[1],reses1[-i-2].name[1],reses1[-i-1].name[1]]
        #check point: all bps are canonical bps
        #make pdb file: center puring in the 1st strand
        skip = False
        for j in range(5):
            if (reses0[i+j].name == "DA" and reses1[-i-j-1].name == "DT") or (reses0[i+j].name == "DT" and reses1[-i-j-1].name == "DA") or (reses0[i+j].name == "DG" and reses1[-i-j-1].name == "DC") or (reses0[i+j].name == "DC" and reses1[-i-j-1].name == "DG"):
                if swap == False:
                    newseg0.reses[j].atoms = reses0[i+j].atoms
                    newseg0.reses[j].name = reses0[i+j].name
                    newseg0.reses[j].resi = j+1
                    newseg1.reses[j].atoms = reses1[-i+j-5].atoms
                    newseg1.reses[j].name = reses1[-i+j-5].name
                    newseg1.reses[j].resi = j+6
                elif swap == True: 
                    newseg0.reses[j].atoms = reses1[-i+j-5].atoms
                    newseg0.reses[j].name = reses1[-i+j-5].name
                    newseg0.reses[j].resi = j+1
                    newseg1.reses[j].atoms = reses0[i+j].atoms
                    newseg1.reses[j].name = reses0[i+j].name
                    newseg1.reses[j].resi = j+6
            else:
                print("=================================")
                print("Unexpected base pair motif! Skip!")
                print("=================================")
                print("Position: bp0 "+str(reses0[i+j].resi)+" "+reses0[i+j].name+" bp1 "+str(reses1[-i-j-1].resi)+" "+reses1[-i-j-1].name)
                skip = True
                break
        if skip: continue
        #parse json information
        if swap == False:
            hform1 = najson.json_file['stems'][stemid_int]['helix_form'][i+1]
            mingw1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['groove_widths'][0]
            majgw1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['groove_widths'][2]
            slide1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][1]
            roll1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][4]
            twist1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][5]
            incli1 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['heli_params'][3]
            hform2 = najson.json_file['stems'][stemid_int]['helix_form'][i+2]
            mingw2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['groove_widths'][0]
            majgw2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['groove_widths'][2]
            slide2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][1]
            roll2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][4]
            twist2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][5]
            incli2 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['heli_params'][3]
        elif swap == True:
            hform1 = najson.json_file['stems'][stemid_int]['helix_form'][i+2]
            mingw1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['groove_widths'][0]
            majgw1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['groove_widths'][2]
            slide1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][1]
            roll1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][4]
            twist1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['step_params'][5]
            incli1 = najson.json_file['stems'][stemid_int]['pairs'][i+2]['heli_params'][3]
            hform2 = najson.json_file['stems'][stemid_int]['helix_form'][i+1]
            mingw2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['groove_widths'][0]
            majgw2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['groove_widths'][2]
            slide2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][1]
            roll2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][4]
            twist2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['step_params'][5]
            incli2 = najson.json_file['stems'][stemid_int]['pairs'][i+1]['heli_params'][3]
        newseg0.chid = "A"
        newseg1.chid = "B"
        newmol.write(oupdir+"/"+unitid+".pdb")
        ouplist.append([unitid,pdbid,stemid,sform,cbpmf,''.join(seq[1:4]),''.join(seq),hform1,hform2,mingw1,mingw2,majgw1,majgw2,slide1,slide2,roll1,roll2,twist1,twist2,incli1,incli2])
        
ouplist.sort()
file = open('Output_seq.txt','wt')
file.write("unitid pdbid stemid sf motif 3bps 5bps hf1 hf2 mingw1 mingw2 majgw1 majgw2 slide1 slide2 roll1 roll2 twist1 twist2 incli1 incli2\n")
for oup in ouplist:
    file.write("%30s %4s %2s %1s %2s %3s %5s %1s %1s %6s %6s %6s %6s %7s %7s %7s %7s %7s %7s %7s %7s\n"%(oup[0],oup[1],oup[2],oup[3],oup[4],oup[5],oup[6],oup[7],oup[8],oup[9],oup[10],oup[11],oup[12],oup[13],oup[14],oup[15],oup[16],oup[17],oup[18],oup[19],oup[20]))

file.close()
toc = time.time()  #timer: end
print(">>>>>> Total time used %f >>>>>>"%(toc-tic))
