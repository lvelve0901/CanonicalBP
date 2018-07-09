#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
from commontool import read, plot2D, hist2D

abg_init = pd.read_table('Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('Output_seq.txt',delim_whitespace=True)
nts_init = pd.read_table('Output_nts.txt',delim_whitespace=True)
nts3_init = pd.read_table('Output_nts_pur+1.txt',delim_whitespace=True)
nts5_init = pd.read_table('Output_nts_pur-1.txt',delim_whitespace=True)
sug_init = pd.read_table('Output_sugar.txt',delim_whitespace=True)


abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]
nts = nts_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]
nts3 = nts3_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]
nts5 = nts5_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]
sug = sug_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B')]

print(len(abg))

fig = pl.figure(1,figsize=(24,7))

hist2D(fig,3,8,1,nts.pur_e_z,abg.beta,10,2,0.002,"","pur_e_z","betah",-180,180,0,60)
hist2D(fig,3,8,2,nts.pur_phase,abg.beta,10,2,0.002,"","pur_phase","betah",0,360,0,60)
hist2D(fig,3,8,3,nts.pur_gamma,abg.beta,10,2,0.002,"","pur_gamma","betah",-180,180,0,60)
hist2D(fig,3,8,4,nts5.pur_e_z,abg.beta,10,2,0.002,"","5'pur_e_z","betah",-180,180,0,60)
hist2D(fig,3,8,5,nts.pym_e_z,abg.beta,10,2,0.002,"","pym_e_z","betah",-180,180,0,60)
hist2D(fig,3,8,6,nts.pym_phase,abg.beta,10,2,0.002,"","pym_phase","betah",0,360,0,60)
hist2D(fig,3,8,7,nts.pym_gamma,abg.beta,10,2,0.002,"","pym_gamma","betah",-180,180,0,60)
hist2D(fig,3,8,8,nts5.pym_e_z,abg.beta,10,2,0.002,"","5'pym_e_z","betah",-180,180,0,60)
#hist2D(fig,3,8,8,sug.bp_c1pc1p,abg.beta,0.5,2,0.002,"","5'pym_e_z","betah",7,14,0,60)


hist2D(fig,3,8,9,nts.pur_e_z,abg.gamma,10,10,0.0002,"","pur_e_z","gammah",-180,180,-180,180)
hist2D(fig,3,8,10,nts.pur_phase,abg.gamma,10,10,0.0002,"","pur_phase","gammah",0,360,-180,180)
hist2D(fig,3,8,11,nts.pur_gamma,abg.gamma,10,10,0.0002,"","pur_gamma","gammah",-180,180,-180,180)
hist2D(fig,3,8,12,nts5.pur_e_z,abg.gamma,10,10,0.0002,"","5'pur_e_z","gammah",-180,180,-180,180)
hist2D(fig,3,8,13,nts.pym_e_z,abg.gamma,10,10,0.0002,"","pym_e_z","gammah",-180,180,-180,180)
hist2D(fig,3,8,14,nts.pym_phase,abg.gamma,10,10,0.0002,"","pym_phase","gammah",0,360,-180,180)
hist2D(fig,3,8,15,nts.pym_gamma,abg.gamma,10,10,0.0002,"","pym_gamma","gammah",-180,180,-180,180)
hist2D(fig,3,8,16,nts5.pym_e_z,abg.gamma,10,10,0.0002,"","5'pym_e_z","gammah",-180,180,-180,180)
#hist2D(fig,3,8,16,sug.bp_c1pc1p,abg.gamma,0.5,10,0.002,"","5'pym_e_z","gammah",7,14,-180,180)

hist2D(fig,3,8,17,nts.pur_e_z,abg.zeta,10,2,0.002,"","pur_e_z","zetah",-180,180,0,60)
hist2D(fig,3,8,18,nts.pur_phase,abg.zeta,10,2,0.002,"","pur_phase","zetah",0,360,0,60)
hist2D(fig,3,8,19,nts.pur_gamma,abg.zeta,10,2,0.002,"","pur_gamma","zetah",-180,180,0,60)
hist2D(fig,3,8,20,nts5.pur_e_z,abg.zeta,10,2,0.002,"","5'pur_e_z","zetah",-180,180,0,60)
hist2D(fig,3,8,21,nts.pym_e_z,abg.zeta,10,2,0.002,"","pym_e_z","zetah",-180,180,0,60)
hist2D(fig,3,8,22,nts.pym_phase,abg.zeta,10,2,0.002,"","pym_phase","zetah",0,360,0,60)
hist2D(fig,3,8,23,nts.pym_gamma,abg.zeta,10,2,0.002,"","pym_gamma","zetah",-180,180,0,60)
hist2D(fig,3,8,24,nts5.pym_e_z,abg.zeta,10,2,0.002,"","5'pym_e_z","zetah",-180,180,0,60)
#hist2D(fig,3,8,24,sug.bp_c1pc1p,abg.zeta,0.5,2,0.002,"","5'pym_e_z","zetah",7,14,0,60)


fig.tight_layout()

pl.show()
