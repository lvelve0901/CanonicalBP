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

fig = pl.figure(1,figsize=(17,14))

hist2D(fig,6,7,1 ,nts.pur_e_z,abg.beta,10,2,0.001,"","pur_e_z","betah",-180,180,0,60)
hist2D(fig,6,7,2 ,nts.pur_epsilon,abg.beta,10,2,0.001,"","pur_epsilon","",-180,180,0,60)
hist2D(fig,6,7,3 ,nts.pur_phase,abg.beta,10,2,0.001,"","pur_phase","",0,360,0,60)
hist2D(fig,6,7,4 ,nts3.pur_a_g,abg.beta,10,2,0.001,"","pur_a_g","",-180,180,0,60)
hist2D(fig,6,7,5 ,nts.pur_gamma,abg.beta,10,2,0.001,"","pur_gamma","",-180,180,0,60)
hist2D(fig,6,7,6 ,nts5.pur_e_z,abg.beta,10,2,0.001,"","5'pur_e_z","",-180,180,0,60)
hist2D(fig,6,7,7 ,nts3.pur_a_g,abg.beta,10,2,0.001,"","3'pur_a_g","",-180,180,0,60)

hist2D(fig,6,7,8 ,nts.pym_e_z,abg.beta,10,2,0.001,"","pym_e_z","betah",-180,180,0,60)
hist2D(fig,6,7,9 ,nts.pym_epsilon,abg.beta,10,2,0.001,"","pym_epsilon","",-180,180,0,60)
hist2D(fig,6,7,10,nts.pym_phase,abg.beta,10,2,0.001,"","pym_phase","",0,360,0,60)
hist2D(fig,6,7,11,nts3.pym_a_g,abg.beta,10,2,0.001,"","pym_a_g","",-180,180,0,60)
hist2D(fig,6,7,12,nts.pym_gamma,abg.beta,10,2,0.001,"","pym_gamma","",-180,180,0,60)
hist2D(fig,6,7,13,nts5.pym_e_z,abg.beta,10,2,0.001,"","5'pym_e_z","",-180,180,0,60)
hist2D(fig,6,7,14,nts3.pym_a_g,abg.beta,10,2,0.001,"","3'pym_a_g","",-180,180,0,60)

hist2D(fig,6,7,15,nts.pur_e_z,abg.gamma,10,10,0.0004,"","pur_e_z","gammah",-180,180,-180,180)
hist2D(fig,6,7,16,nts.pur_epsilon,abg.gamma,10,10,0.0004,"","pur_epsilon","",-180,180,-180,180)
hist2D(fig,6,7,17,nts.pur_phase,abg.gamma,10,10,0.0004,"","pur_phase","",0,360,-180,180)
hist2D(fig,6,7,18,nts3.pur_a_g,abg.gamma,10,10,0.0004,"","pur_a_g","",-180,180,-180,180)
hist2D(fig,6,7,19,nts.pur_gamma,abg.gamma,10,10,0.0004,"","pur_gamma","",-180,180,-180,180)
hist2D(fig,6,7,20,nts5.pur_e_z,abg.gamma,10,10,0.0004,"","5'pur_e_z","",-180,180,-180,180)
hist2D(fig,6,7,21,nts3.pur_a_g,abg.gamma,10,10,0.0004,"","3'pur_a_g","",-180,180,-180,180)

hist2D(fig,6,7,22,nts.pym_e_z,abg.gamma,10,10,0.0004,"","pym_e_z","gammah",-180,180,-180,180)
hist2D(fig,6,7,23,nts.pym_epsilon,abg.gamma,10,10,0.0004,"","pym_epsilon","",-180,180,-180,180)
hist2D(fig,6,7,24,nts.pym_phase,abg.gamma,10,10,0.0004,"","pym_phase","",0,360,-180,180)
hist2D(fig,6,7,25,nts3.pym_a_g,abg.gamma,10,10,0.0004,"","pym_a_g","",-180,180,-180,180)
hist2D(fig,6,7,26,nts.pym_gamma,abg.gamma,10,10,0.0004,"","pym_gamma","",-180,180,-180,180)
hist2D(fig,6,7,27,nts5.pym_e_z,abg.gamma,10,10,0.0004,"","5'pym_e_z","",-180,180,-180,180)
hist2D(fig,6,7,28,nts3.pym_a_g,abg.gamma,10,10,0.0004,"","3'pym_a_g","",-180,180,-180,180)

hist2D(fig,6,7,29,nts.pur_e_z,abg.zeta,10,2,0.001,"","pur_e_z","zetah",-180,180,0,60)
hist2D(fig,6,7,30,nts.pur_epsilon,abg.zeta,10,2,0.001,"","pur_epsilon","",-180,180,0,60)
hist2D(fig,6,7,31,nts.pur_phase,abg.zeta,10,2,0.001,"","pur_phase","",0,360,0,60)
hist2D(fig,6,7,32,nts3.pur_a_g,abg.zeta,10,2,0.001,"","pur_a_g","",-180,180,0,60)
hist2D(fig,6,7,33,nts.pur_gamma,abg.zeta,10,2,0.001,"","pur_gamma","",-180,180,0,60)
hist2D(fig,6,7,34,nts5.pur_e_z,abg.zeta,10,2,0.001,"","5'pur_e_z","",-180,180,0,60)
hist2D(fig,6,7,35,nts3.pur_a_g,abg.zeta,10,2,0.001,"","3'pur_a_g","",-180,180,0,60)

hist2D(fig,6,7,36,nts.pym_e_z,abg.zeta,10,2,0.001,"","pym_e_z","zetah",-180,180,0,60)
hist2D(fig,6,7,37,nts.pym_epsilon,abg.zeta,10,2,0.001,"","pym_epsilon","",-180,180,0,60)
hist2D(fig,6,7,38,nts.pym_phase,abg.zeta,10,2,0.001,"","pym_phase","",0,360,0,60)
hist2D(fig,6,7,39,nts3.pym_a_g,abg.zeta,10,2,0.001,"","pym_a_g","",-180,180,0,60)
hist2D(fig,6,7,40,nts.pym_gamma,abg.zeta,10,2,0.001,"","pym_gamma","",-180,180,0,60)
hist2D(fig,6,7,41,nts5.pym_e_z,abg.zeta,10,2,0.001,"","5'pym_e_z","",-180,180,0,60)
hist2D(fig,6,7,42,nts3.pym_a_g,abg.zeta,10,2,0.001,"","3'pym_a_g","",-180,180,0,60)


fig.tight_layout()

fig = pl.figure(2,figsize=(21,7))

hist2D(fig,3,7,1 ,seq.mingw1.convert_objects(convert_numeric=True),abg.beta,0.5,2,0.01,"","minor_gw1","betah",0,20,0,60)
hist2D(fig,3,7,2 ,seq.mingw2.convert_objects(convert_numeric=True),abg.beta,0.5,2,0.01,"","minor_gw2","betah",0,20,0,60)
hist2D(fig,3,7,3 ,seq.roll1.convert_objects(convert_numeric=True),abg.beta,2,2,0.004,"","roll1","betah",-40,40,0,60)
hist2D(fig,3,7,4 ,seq.roll1.convert_objects(convert_numeric=True),abg.beta,2,2,0.004,"","roll2","betah",-40,40,0,60)
hist2D(fig,3,7,5 ,seq.twist1.convert_objects(convert_numeric=True),abg.beta,2,2,0.004,"","twist1","betah",10,70,0,60)
hist2D(fig,3,7,6 ,seq.twist2.convert_objects(convert_numeric=True),abg.beta,2,2,0.004,"","twist2","betah",10,70,0,60)
hist2D(fig,3,7,7 ,sug.bp_c1pc1p,abg.beta,0.2,2,0.02,"","c1pc1p","betah",6,14,0,60)

hist2D(fig,3,7,8 ,seq.mingw1.convert_objects(convert_numeric=True),abg.gamma,0.5,10,0.001,"","minor_gw1","gammah",0,20,-180,180)
hist2D(fig,3,7,9 ,seq.mingw2.convert_objects(convert_numeric=True),abg.gamma,0.5,10,0.001,"","minor_gw2","gammah",0,20,-180,180)
hist2D(fig,3,7,10,seq.roll1.convert_objects(convert_numeric=True),abg.gamma,2,10,0.0004,"","roll1","gammah",-40,40,-180,180)
hist2D(fig,3,7,11,seq.roll1.convert_objects(convert_numeric=True),abg.gamma,2,10,0.0004,"","roll2","gammah",-40,40,-180,180)
hist2D(fig,3,7,12,seq.twist1.convert_objects(convert_numeric=True),abg.gamma,2,10,0.0004,"","twist1","gammah",10,70,-180,180)
hist2D(fig,3,7,13,seq.twist2.convert_objects(convert_numeric=True),abg.gamma,2,10,0.0004,"","twist2","gammah",10,70,-180,180)
hist2D(fig,3,7,14,sug.bp_c1pc1p,abg.gamma,0.2,10,0.002,"","c1pc1p","gammah",6,14,-180,180)

hist2D(fig,3,7,15,seq.mingw1.convert_objects(convert_numeric=True),abg.zeta,0.5,2,0.01,"","minor_gw1","zetah",0,20,0,60)
hist2D(fig,3,7,16,seq.mingw2.convert_objects(convert_numeric=True),abg.zeta,0.5,2,0.01,"","minor_gw2","zetah",0,20,0,60)
hist2D(fig,3,7,17,seq.roll1.convert_objects(convert_numeric=True),abg.zeta,2,2,0.004,"","roll1","zetah",-40,40,0,60)
hist2D(fig,3,7,18,seq.roll1.convert_objects(convert_numeric=True),abg.zeta,2,2,0.004,"","roll2","zetah",-40,40,0,60)
hist2D(fig,3,7,19,seq.twist1.convert_objects(convert_numeric=True),abg.zeta,2,2,0.004,"","twist1","zetah",10,70,0,60)
hist2D(fig,3,7,20,seq.twist2.convert_objects(convert_numeric=True),abg.zeta,2,2,0.004,"","twist2","zetah",10,70,0,60)
hist2D(fig,3,7,21,sug.bp_c1pc1p,abg.zeta,0.2,2,0.02,"","c1pc1p","zetah",6,14,0,60)

fig.tight_layout()


pl.show()
