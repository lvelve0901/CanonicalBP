#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
from commontool import read, plot2D, hist2D

abg_init = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('../Output_seq.txt',delim_whitespace=True)
nts_init = pd.read_table('../Output_nts.txt',delim_whitespace=True)
nts3_init = pd.read_table('../Output_nts_pur+1.txt',delim_whitespace=True)
nts5_init = pd.read_table('../Output_nts_pur-1.txt',delim_whitespace=True)
sug_init = pd.read_table('../Output_sugar.txt',delim_whitespace=True)
l = [i for s in read('dna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
nts = nts_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
nts3 = nts3_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
nts5 = nts5_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
sug = sug_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]

print(len(abg))

fig = pl.figure(1,figsize=(24,7))

hist2D(fig,3,8,1 ,nts.pur_e_z,nts.pur_phase,10,10,0.0001,"","pur_e_z","pur_phase",-180,180,0,360)
hist2D(fig,3,8,2 ,nts.pur_e_z,nts.pur_gamma,10,10,0.0001,"","pur_e_z","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,3 ,nts.pur_e_z,nts.pur_a_g,10,10,0.0001,"","pur_e_z","pur_a_g",-180,180,-180,180)
hist2D(fig,3,8,4 ,nts.pur_epsilon,nts.pur_phase,10,10,0.0001,"","pur_epsilon","pur_phase",-180,180,0,360)
hist2D(fig,3,8,5 ,nts.pur_epsilon,nts.pur_gamma,10,10,0.0001,"","pur_epsilon","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,6 ,nts.pur_gamma,nts.pur_phase,10,10,0.0001,"","pur_gamma","pur_phase",-180,180,0,360)
hist2D(fig,3,8,7 ,nts.pur_chi,nts.pur_phase,10,10,0.0001,"","pur_chi","pur_phase",-180,180,0,360)
hist2D(fig,3,8,8 ,nts.pur_alpha,nts.pur_gamma,10,10,0.0001,"","pur_alpha","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,9 ,nts.pur_epsilon,nts.pur_zeta,10,10,0.0001,"","pur_epsilon","pur_zeta",-180,180,-180,180)
hist2D(fig,3,8,10,nts.pur_e_z,nts.pur_epsilon,10,10,0.0001,"","pur_e_z","pur_epsilon",-180,180,-180,180)
hist2D(fig,3,8,11,nts.pur_e_z,nts.pur_zeta,10,10,0.0001,"","pur_e_z","pur_zeta",-180,180,-180,180)
hist2D(fig,3,8,12,nts.pur_a_g,nts.pur_alpha,10,10,0.0001,"","pur_a_g","pur_alpha",-180,180,-180,180)
hist2D(fig,3,8,13,nts.pur_a_g,nts.pur_gamma,10,10,0.0001,"","pur_a_g","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,14,nts.pur_e_z,nts.pym_e_z,10,10,0.0001,"","pur_e_z","pym_e_z",-180,180,-180,180)
hist2D(fig,3,8,15,nts.pur_a_g,nts.pym_a_g,10,10,0.0001,"","pur_a_g","pym_a_g",-180,180,-180,180)
hist2D(fig,3,8,16,nts.pur_e_z,seq.mingw1.convert_objects(convert_numeric=True),10,0.5,0.001,"","pur_e_z","mingw1",-180,180,5,25)
hist2D(fig,3,8,17,nts.pur_a_g,seq.mingw1.convert_objects(convert_numeric=True),10,0.5,0.001,"","pur_a_g","mingw1",-180,180,5,25)
hist2D(fig,3,8,18,nts.pur_phase,seq.mingw1.convert_objects(convert_numeric=True),10,0.5,0.001,"","pur_phase","mingw1",0,360,5,25)
hist2D(fig,3,8,19,seq.slide1.convert_objects(convert_numeric=True),seq.roll1.convert_objects(convert_numeric=True),0.2,1,0.06,"","slide1","roll1",-2,4,-30,30)
hist2D(fig,3,8,20,nts.pur_gamma,seq.mingw1.convert_objects(convert_numeric=True),10,0.5,0.001,"","pur_gamma","mingw1",-180,180,5,25)
hist2D(fig,3,8,21,sug.bp_c1pc1p,seq.twist1.convert_objects(convert_numeric=True),0.1,1,0.1,"","bp_c1pc1p","twist",8,14,10,50)
hist2D(fig,3,8,22,nts.pur_chi,seq.twist1,3,1,0.003,"","pur_chi","twist1",-180,0,0,60)


fig.tight_layout()

pl.show()
