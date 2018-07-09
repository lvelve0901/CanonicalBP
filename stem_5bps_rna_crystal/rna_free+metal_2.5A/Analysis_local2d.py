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
l = [i for s in read('rna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
nts = nts_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
nts3 = nts3_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
nts5 = nts5_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
sug = sug_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]

print(len(abg))

fig = pl.figure(1,figsize=(24,7))

hist2D(fig,3,8,1 ,nts.pur_e_z,nts.pur_phase,10,10,0.0001,"","pur_e_z","pur_phase",-180,180,0,360)
hist2D(fig,3,8,2 ,nts.pur_e_z,nts.pur_gamma,10,10,0.0001,"","pur_e_z","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,3 ,nts.pur_e_z,nts.pur_a_g,10,10,0.0001,"","pur_e_z","pur_a_g",-180,180,-180,180)
hist2D(fig,3,8,4 ,nts.pur_epsilon,nts.pur_phase,10,10,0.0001,"","pur_epsilon","pur_phase",-180,180,0,360)
hist2D(fig,3,8,5 ,nts.pur_epsilon,nts.pur_gamma,10,10,0.0001,"","pur_epsilon","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,6 ,nts.pur_gamma,nts.pur_phase,10,10,0.0001,"","pur_gamma","pur_phase",-180,180,0,360)
hist2D(fig,3,8,7 ,nts.pur_chi,nts.pur_phase,1,1,0.01,"","pur_chi","pur_phase",-180,-90,0,60)
hist2D(fig,3,8,8 ,nts.pur_alpha,nts.pur_gamma,10,10,0.0001,"","pur_alpha","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,9 ,nts.pur_epsilon,nts.pur_zeta,10,10,0.0001,"","pur_epsilon","pur_zeta",-180,180,-180,180)
hist2D(fig,3,8,10,nts.pur_e_z,nts.pur_epsilon,10,10,0.0001,"","pur_e_z","pur_epsilon",-180,180,-180,180)
hist2D(fig,3,8,11,nts.pur_e_z,nts.pur_zeta,10,10,0.0001,"","pur_e_z","pur_zeta",-180,180,-180,180)
hist2D(fig,3,8,12,nts.pur_a_g,nts.pur_alpha,10,10,0.0001,"","pur_a_g","pur_alpha",-180,180,-180,180)
hist2D(fig,3,8,13,nts.pur_a_g,nts.pur_gamma,10,10,0.0001,"","pur_a_g","pur_gamma",-180,180,-180,180)
hist2D(fig,3,8,14,nts.pur_chi,seq.roll1,1,2,0.01,"","pur_chi","roll1",-180,-90,-30,30)
hist2D(fig,3,8,15,nts.pur_chi,seq.twist1,1,1,0.01,"","pur_chi","twist1",-180,-90,0,60)


fig.tight_layout()

pl.show()
