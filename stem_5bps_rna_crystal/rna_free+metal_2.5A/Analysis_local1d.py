#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
from commontool import read, plot2D, hist1D, hist2D

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

fig = pl.figure(1,figsize=(32,24))

hist1D(fig,4,8,1 , nts.pur_alpha,10,'alpha','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,2 , nts.pur_beta,10,'beta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,3 , nts.pur_gamma,10,'gamma','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,4 , nts.pur_delta,10,'delta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,5 , nts.pur_epsilon,10,'epsilon','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,6 , nts.pur_zeta,10,'zeta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,7 , nts.pur_a_g,10,'a-g','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,8 , nts.pur_e_z,10,'e-z','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,9 , nts.pur_chi,10,'chi','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,10, nts.pur_phase,10,'phase','red',0,360,0.2,legend=None)
hist1D(fig,4,8,11, nts.pur_ampli,2,'ampli','red',0,60,0.2,legend=None)
hist1D(fig,4,8,12, nts.pur_v0,10,'v0','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,13, nts.pur_v1,10,'v1','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,14, nts.pur_v2,10,'v2','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,15, nts.pur_v3,10,'v3','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,16, nts.pur_v4,10,'v4','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,17, nts.pym_alpha,10,'alpha','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,18, nts.pym_beta,10,'beta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,19, nts.pym_gamma,10,'gamma','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,20, nts.pym_delta,10,'delta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,21, nts.pym_epsilon,10,'epsilon','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,22, nts.pym_zeta,10,'zeta','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,23, nts.pym_a_g,10,'a-g','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,24, nts.pym_e_z,10,'e-z','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,25, nts.pym_chi,10,'chi','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,26, nts.pym_phase,10,'phase','red',0,360,0.2,legend=None)
hist1D(fig,4,8,27, nts.pym_ampli,2,'ampli','red',0,60,0.2,legend=None)
hist1D(fig,4,8,28, nts.pym_v0,10,'v0','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,29, nts.pym_v1,10,'v1','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,30, nts.pym_v2,10,'v2','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,31, nts.pym_v3,10,'v3','red',-180,180,0.2,legend=None)
hist1D(fig,4,8,32, nts.pym_v4,10,'v4','red',-180,180,0.2,legend=None)

fig.tight_layout()

fig = pl.figure(2,figsize=(32,24))

hist1D(fig,4,8,1 , seq.roll1,2,'roll1','red',-30,30,0.2,legend=None)
hist1D(fig,4,8,2 , seq.twist1,2,'twist1','red',0,60,0.2,legend=None)
hist1D(fig,4,8,3 , seq.incli1,2,'incli1','red',-30,30,0.2,legend=None)
x = seq.mingw1.convert_objects(convert_numeric=True)
hist1D(fig,4,8,4 , x[~np.isnan(x)],0.5,'mingw1','red',5,25,0.2,legend=None)

hist1D(fig,4,8,9 , seq.roll2,2,'roll2','red',-30,30,0.2,legend=None)
hist1D(fig,4,8,10, seq.twist2,2,'twist2','red',0,60,0.2,legend=None)
hist1D(fig,4,8,11, seq.incli2,2,'incli2','red',-30,30,0.2,legend=None)
x = seq.mingw2.convert_objects(convert_numeric=True)
hist1D(fig,4,8,12, x[~np.isnan(x)],0.5,'mingw2','red',5,25,0.2,legend=None)

hist1D(fig,4,8,17, sug.bp_c1pc1p,0.5,'c1pc1p','red',5,15,0.2,legend=None)

fig.tight_layout()

pl.show()
