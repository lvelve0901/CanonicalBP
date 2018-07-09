#!/usr/bin/python

import json # Handle JSON DSSR files
import os
import sys
import pandas as pd
import numpy as np
import re
import learnna_json as lna_json
from pdblib.base import *

os.system("mkdir -p Json")
json_dir = './Json'
pdb_dir = 'Fragment_nts'
oupname = 'Output_nts_pur-1'

#change range of -360~-180 --> 0~180 and 180~360 --> -180~0
def restr(number):
    if number > 180:
        number = number - 360
    elif number < -180:
        number = number + 360
    return number

pdblist = []  #list to store PDB file names
for file in os.listdir(pdb_dir): #read input directory
    if file.endswith(".pdb"):
        pdblist.append(file)

pdblist.sort()
output_nts = []  #list to store output information

#Loop over list of files in data dir, process with DSSR to get json file
for idx, pdb in enumerate(pdblist):
    unitid = str(pdb)[:-4]
    pdb_f = os.path.join(pdb_dir,pdb)  #path to pdb file
    json_f = os.path.join(json_dir,pdb.replace(".pdb",".json"))  #path to json file
    print ("--- Working on [%s] (%d of %d) ---"%(pdb,idx+1,len(pdblist)))
    if os.path.isfile(pdb_f):  #DSSR convert PDB to json files
    	os.system("x3dna-dssr --json --more -i=%s -o=%s"%(pdb_f,json_f))
    	os.system("x3dna-dssr --cleanup")
    najson = lna_json.NA_JSON()  #initialize class objects
    with open(json_f) as json_data:  #read each json file
    	data = json.load(json_data)
    najson.set_json(data)  #pass json file to class pbject
    najson.read_idx()  #set index from own json file
    purine = najson.json_file['nts'][1]
    pyrimi = najson.json_file['nts'][8]
    purine_alpha = purine['alpha'] 
    purine_beta = purine['beta'] 
    purine_gamma = purine['gamma'] 
    purine_delta = purine['delta'] 
    purine_epsilon = purine['epsilon'] 
    purine_zeta = purine['zeta'] 
    purine_e_z = restr(purine['epsilon_zeta']) 
    purine_a_g = restr(purine['alpha'] - purine['gamma'])
    purine_chi = purine['chi']
    purine_v0 = purine['v0'] 
    purine_v1 = purine['v1'] 
    purine_v2 = purine['v2'] 
    purine_v3 = purine['v3'] 
    purine_v4 = purine['v4'] 
    purine_phase = purine['phase_angle'] 
    purine_ampli = purine['amplitude'] 
    purine_pucker = purine['puckering'] 
    purine_bb = purine['bb_type'] 
    pyrimi_alpha = pyrimi['alpha'] 
    pyrimi_beta = pyrimi['beta'] 
    pyrimi_gamma = pyrimi['gamma'] 
    pyrimi_delta = pyrimi['delta'] 
    pyrimi_epsilon = pyrimi['epsilon'] 
    pyrimi_zeta = pyrimi['zeta'] 
    pyrimi_e_z = restr(pyrimi['epsilon_zeta']) 
    pyrimi_a_g = restr(pyrimi['alpha'] - pyrimi['gamma'])
    pyrimi_chi = pyrimi['chi']
    pyrimi_v0 = pyrimi['v0'] 
    pyrimi_v1 = pyrimi['v1'] 
    pyrimi_v2 = pyrimi['v2'] 
    pyrimi_v3 = pyrimi['v3'] 
    pyrimi_v4 = pyrimi['v4'] 
    pyrimi_phase = pyrimi['phase_angle'] 
    pyrimi_ampli = pyrimi['amplitude'] 
    pyrimi_pucker = pyrimi['puckering'] 
    pyrimi_bb = pyrimi['bb_type'] 
    output_nts.append([unitid,purine_alpha,purine_beta,purine_gamma,purine_delta,purine_epsilon,purine_zeta,purine_e_z,purine_a_g,purine_chi,purine_v0,purine_v1,purine_v2,purine_v3,purine_v4,purine_phase,purine_ampli,purine_pucker,purine_bb,pyrimi_alpha,pyrimi_beta,pyrimi_gamma,pyrimi_delta,pyrimi_epsilon,pyrimi_zeta,pyrimi_e_z,pyrimi_a_g,pyrimi_chi,pyrimi_v0,pyrimi_v1,pyrimi_v2,pyrimi_v3,pyrimi_v4,pyrimi_phase,pyrimi_ampli,pyrimi_pucker,pyrimi_bb])

file_nts = open('%s.txt'%oupname,'w')

file_nts.write('unitid pur_alpha pur_beta pur_gamma pur_delta pur_epsilon pur_zeta pur_e_z pur_a_g pur_chi pur_v0 pur_v1 pur_v2 pur_v3 pur_v4 pur_phase pur_ampli pur_pucker pur_bb pym_alpha pym_beta pym_gamma pym_delta pym_epsilon pym_zeta pym_e_z pym_a_g pym_chi pym_v0 pym_v1 pym_v2 pym_v3 pym_v4 pym_phase pym_ampli pym_pucker pym_bb\n')
for nt in output_nts:
	file_nts.write('%30s %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %9s %3s %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %9s %3s\n'%(nt[0],nt[1],nt[2],nt[3],nt[4],nt[5],nt[6],nt[7],nt[8],nt[9],nt[10],nt[11],nt[12],nt[13],nt[14],nt[15],nt[16],nt[17],nt[18],nt[19],nt[20],nt[21],nt[22],nt[23],nt[24],nt[25],nt[26],nt[27],nt[28],nt[29],nt[30],nt[31],nt[32],nt[33],nt[34],nt[35],nt[36]))

file_nts.close()
