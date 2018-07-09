
#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
import matplotlib.backends.backend_pdf
from commontool import read, plot2D, hist1D, hist2D, phist2D

name = "Plotseq_dna_free+metal_XPX.pdf"

abg_init = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('../Output_seq.txt',delim_whitespace=True)
l = [i for s in read('dna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]

AAA = abg.loc[seq['3bps'] == 'AAA']
AAG = abg.loc[seq['3bps'] == 'AAG']
AAC = abg.loc[seq['3bps'] == 'AAC']
AAT = abg.loc[seq['3bps'] == 'AAT']
GAA = abg.loc[seq['3bps'] == 'GAA']
GAG = abg.loc[seq['3bps'] == 'GAG']
GAC = abg.loc[seq['3bps'] == 'GAC']
GAT = abg.loc[seq['3bps'] == 'GAT']
CAA = abg.loc[seq['3bps'] == 'CAA']
CAG = abg.loc[seq['3bps'] == 'CAG']
CAC = abg.loc[seq['3bps'] == 'CAC']
CAT = abg.loc[seq['3bps'] == 'CAT']
TAA = abg.loc[seq['3bps'] == 'TAA']
TAG = abg.loc[seq['3bps'] == 'TAG']
TAC = abg.loc[seq['3bps'] == 'TAC']
TAT = abg.loc[seq['3bps'] == 'TAT']
AGA = abg.loc[seq['3bps'] == 'AGA']
AGG = abg.loc[seq['3bps'] == 'AGG']
AGC = abg.loc[seq['3bps'] == 'AGC']
AGT = abg.loc[seq['3bps'] == 'AGT']
GGA = abg.loc[seq['3bps'] == 'GGA']
GGG = abg.loc[seq['3bps'] == 'GGG']
GGC = abg.loc[seq['3bps'] == 'GGC']
GGT = abg.loc[seq['3bps'] == 'GGT']
CGA = abg.loc[seq['3bps'] == 'CGA']
CGG = abg.loc[seq['3bps'] == 'CGG']
CGC = abg.loc[seq['3bps'] == 'CGC']
CGT = abg.loc[seq['3bps'] == 'CGT']
TGA = abg.loc[seq['3bps'] == 'TGA']
TGG = abg.loc[seq['3bps'] == 'TGG']
TGC = abg.loc[seq['3bps'] == 'TGC']
TGT = abg.loc[seq['3bps'] == 'TGT']

fig1 = pl.figure(1,figsize=(27,13))

fig1.suptitle('beta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.beta,1,'AAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,2 , AAG.beta,1,'AAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,3 , AAC.beta,1,'AAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,4 , AAT.beta,1,'AAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,5 , GAA.beta,1,'GAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,6 , GAG.beta,1,'GAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,7 , GAC.beta,1,'GAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,8 , GAT.beta,1,'GAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,9 , CAA.beta,1,'CAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,10, CAG.beta,1,'CAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,11, CAC.beta,1,'CAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,12, CAT.beta,1,'CAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,13, TAA.beta,1,'TAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,14, TAG.beta,1,'TAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,15, TAC.beta,1,'TAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,16, TAT.beta,1,'TAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,17, AGA.beta,1,'AGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,18, AGG.beta,1,'AGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,19, AGC.beta,1,'AGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,20, AGT.beta,1,'AGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,21, GGA.beta,1,'GGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,22, GGG.beta,1,'GGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,23, GGC.beta,1,'GGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,24, GGT.beta,1,'GGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,25, CGA.beta,1,'CGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,26, CGG.beta,1,'CGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,27, CGC.beta,1,'CGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,28, CGT.beta,1,'CGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,29, TGA.beta,1,'TGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,30, TGG.beta,1,'TGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,31, TGC.beta,1,'TGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,32, TGT.beta,1,'TGT','red',0,60,0.2,legend=None)

fig1 = pl.figure(2,figsize=(27,13))

fig1.suptitle('gamma distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.gamma,6,'AAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,2 , AAG.gamma,6,'AAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,3 , AAC.gamma,6,'AAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,4 , AAT.gamma,6,'AAT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,5 , GAA.gamma,6,'GAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,6 , GAG.gamma,6,'GAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,7 , GAC.gamma,6,'GAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,8 , GAT.gamma,6,'GAT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,9 , CAA.gamma,6,'CAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,10, CAG.gamma,6,'CAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,11, CAC.gamma,6,'CAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,12, CAT.gamma,6,'CAT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,13, TAA.gamma,6,'TAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,14, TAG.gamma,6,'TAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,15, TAC.gamma,6,'TAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,16, TAT.gamma,6,'TAT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,17, AGA.gamma,6,'AGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,18, AGG.gamma,6,'AGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,19, AGC.gamma,6,'AGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,20, AGT.gamma,6,'AGT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,21, GGA.gamma,6,'GGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,22, GGG.gamma,6,'GGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,23, GGC.gamma,6,'GGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,24, GGT.gamma,6,'GGT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,25, CGA.gamma,6,'CGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,26, CGG.gamma,6,'CGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,27, CGC.gamma,6,'CGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,28, CGT.gamma,6,'CGT','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,29, TGA.gamma,6,'TGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,30, TGG.gamma,6,'TGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,31, TGC.gamma,6,'TGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,32, TGT.gamma,6,'TGT','red',-180,180,0.1,legend=None)

fig1 = pl.figure(3,figsize=(27,13))

fig1.suptitle('zeta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.zeta,1,'AAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,2 , AAG.zeta,1,'AAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,3 , AAC.zeta,1,'AAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,4 , AAT.zeta,1,'AAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,5 , GAA.zeta,1,'GAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,6 , GAG.zeta,1,'GAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,7 , GAC.zeta,1,'GAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,8 , GAT.zeta,1,'GAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,9 , CAA.zeta,1,'CAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,10, CAG.zeta,1,'CAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,11, CAC.zeta,1,'CAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,12, CAT.zeta,1,'CAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,13, TAA.zeta,1,'TAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,14, TAG.zeta,1,'TAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,15, TAC.zeta,1,'TAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,16, TAT.zeta,1,'TAT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,17, AGA.zeta,1,'AGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,18, AGG.zeta,1,'AGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,19, AGC.zeta,1,'AGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,20, AGT.zeta,1,'AGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,21, GGA.zeta,1,'GGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,22, GGG.zeta,1,'GGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,23, GGC.zeta,1,'GGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,24, GGT.zeta,1,'GGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,25, CGA.zeta,1,'CGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,26, CGG.zeta,1,'CGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,27, CGC.zeta,1,'CGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,28, CGT.zeta,1,'CGT','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,29, TGA.zeta,1,'TGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,30, TGG.zeta,1,'TGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,31, TGC.zeta,1,'TGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,32, TGT.zeta,1,'TGT','red',0,60,0.2,legend=None)

fig1 = pl.figure(4,figsize=(27,13))

fig1.suptitle('gamma-beta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,8,1 ,AAA.gamma*np.pi/180,AAA.beta,np.pi/20,2,0.08,'AAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,2 ,AAG.gamma*np.pi/180,AAG.beta,np.pi/20,2,0.08,'AAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,3 ,AAC.gamma*np.pi/180,AAC.beta,np.pi/20,2,0.08,'AAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,4 ,AAT.gamma*np.pi/180,AAT.beta,np.pi/20,2,0.08,'AAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,5 ,GAA.gamma*np.pi/180,GAA.beta,np.pi/20,2,0.08,'GAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,6 ,GAG.gamma*np.pi/180,GAG.beta,np.pi/20,2,0.08,'GAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,7 ,GAC.gamma*np.pi/180,GAC.beta,np.pi/20,2,0.08,'GAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,8 ,GAT.gamma*np.pi/180,GAT.beta,np.pi/20,2,0.08,'GAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,9 ,CAA.gamma*np.pi/180,CAA.beta,np.pi/20,2,0.08,'CAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,10,CAG.gamma*np.pi/180,CAG.beta,np.pi/20,2,0.08,'CAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,11,CAC.gamma*np.pi/180,CAC.beta,np.pi/20,2,0.08,'CAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,12,CAT.gamma*np.pi/180,CAT.beta,np.pi/20,2,0.08,'CAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,13,TAA.gamma*np.pi/180,TAA.beta,np.pi/20,2,0.08,'TAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,14,TAG.gamma*np.pi/180,TAG.beta,np.pi/20,2,0.08,'TAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,15,TAC.gamma*np.pi/180,TAC.beta,np.pi/20,2,0.08,'TAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,16,TAT.gamma*np.pi/180,TAT.beta,np.pi/20,2,0.08,'TAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,17,AGA.gamma*np.pi/180,AGA.beta,np.pi/20,2,0.08,'AGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,18,AGG.gamma*np.pi/180,AGG.beta,np.pi/20,2,0.08,'AGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,19,AGC.gamma*np.pi/180,AGC.beta,np.pi/20,2,0.08,'AGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,20,AGT.gamma*np.pi/180,AGT.beta,np.pi/20,2,0.08,'AGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,21,GGA.gamma*np.pi/180,GGA.beta,np.pi/20,2,0.08,'GGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,22,GGG.gamma*np.pi/180,GGG.beta,np.pi/20,2,0.08,'GGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,23,GGC.gamma*np.pi/180,GGC.beta,np.pi/20,2,0.08,'GGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,24,GGT.gamma*np.pi/180,GGT.beta,np.pi/20,2,0.08,'GGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,25,CGA.gamma*np.pi/180,CGA.beta,np.pi/20,2,0.08,'CGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,26,CGG.gamma*np.pi/180,CGG.beta,np.pi/20,2,0.08,'CGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,27,CGC.gamma*np.pi/180,CGC.beta,np.pi/20,2,0.08,'CGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,28,CGT.gamma*np.pi/180,CGT.beta,np.pi/20,2,0.08,'CGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,29,TGA.gamma*np.pi/180,TGA.beta,np.pi/20,2,0.08,'TGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,30,TGG.gamma*np.pi/180,TGG.beta,np.pi/20,2,0.08,'TGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,31,TGC.gamma*np.pi/180,TGC.beta,np.pi/20,2,0.08,'TGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,32,TGT.gamma*np.pi/180,TGT.beta,np.pi/20,2,0.08,'TGT',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(5,figsize=(27,13))

fig1.suptitle('gamma-zeta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,8,1 ,AAA.gamma*np.pi/180,AAA.zeta,np.pi/20,2,0.08,'AAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,2 ,AAG.gamma*np.pi/180,AAG.zeta,np.pi/20,2,0.08,'AAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,3 ,AAC.gamma*np.pi/180,AAC.zeta,np.pi/20,2,0.08,'AAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,4 ,AAT.gamma*np.pi/180,AAT.zeta,np.pi/20,2,0.08,'AAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,5 ,GAA.gamma*np.pi/180,GAA.zeta,np.pi/20,2,0.08,'GAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,6 ,GAG.gamma*np.pi/180,GAG.zeta,np.pi/20,2,0.08,'GAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,7 ,GAC.gamma*np.pi/180,GAC.zeta,np.pi/20,2,0.08,'GAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,8 ,GAT.gamma*np.pi/180,GAT.zeta,np.pi/20,2,0.08,'GAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,9 ,CAA.gamma*np.pi/180,CAA.zeta,np.pi/20,2,0.08,'CAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,10,CAG.gamma*np.pi/180,CAG.zeta,np.pi/20,2,0.08,'CAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,11,CAC.gamma*np.pi/180,CAC.zeta,np.pi/20,2,0.08,'CAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,12,CAT.gamma*np.pi/180,CAT.zeta,np.pi/20,2,0.08,'CAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,13,TAA.gamma*np.pi/180,TAA.zeta,np.pi/20,2,0.08,'TAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,14,TAG.gamma*np.pi/180,TAG.zeta,np.pi/20,2,0.08,'TAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,15,TAC.gamma*np.pi/180,TAC.zeta,np.pi/20,2,0.08,'TAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,16,TAT.gamma*np.pi/180,TAT.zeta,np.pi/20,2,0.08,'TAT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,17,AGA.gamma*np.pi/180,AGA.zeta,np.pi/20,2,0.08,'AGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,18,AGG.gamma*np.pi/180,AGG.zeta,np.pi/20,2,0.08,'AGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,19,AGC.gamma*np.pi/180,AGC.zeta,np.pi/20,2,0.08,'AGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,20,AGT.gamma*np.pi/180,AGT.zeta,np.pi/20,2,0.08,'AGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,21,GGA.gamma*np.pi/180,GGA.zeta,np.pi/20,2,0.08,'GGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,22,GGG.gamma*np.pi/180,GGG.zeta,np.pi/20,2,0.08,'GGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,23,GGC.gamma*np.pi/180,GGC.zeta,np.pi/20,2,0.08,'GGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,24,GGT.gamma*np.pi/180,GGT.zeta,np.pi/20,2,0.08,'GGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,25,CGA.gamma*np.pi/180,CGA.zeta,np.pi/20,2,0.08,'CGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,26,CGG.gamma*np.pi/180,CGG.zeta,np.pi/20,2,0.08,'CGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,27,CGC.gamma*np.pi/180,CGC.zeta,np.pi/20,2,0.08,'CGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,28,CGT.gamma*np.pi/180,CGT.zeta,np.pi/20,2,0.08,'CGT',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,29,TGA.gamma*np.pi/180,TGA.zeta,np.pi/20,2,0.08,'TGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,30,TGG.gamma*np.pi/180,TGG.zeta,np.pi/20,2,0.08,'TGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,31,TGC.gamma*np.pi/180,TGC.zeta,np.pi/20,2,0.08,'TGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,32,TGT.gamma*np.pi/180,TGT.zeta,np.pi/20,2,0.08,'TGT',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(6,figsize=(27,13))

fig1.suptitle('beta-zeta distribution',fontsize=14,fontweight='bold')

hist2D(fig1,4,8,1 ,AAA.beta,AAA.zeta,2,2,0.01,'AAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,2 ,AAG.beta,AAG.zeta,2,2,0.01,'AAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,3 ,AAC.beta,AAC.zeta,2,2,0.01,'AAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,4 ,AAT.beta,AAT.zeta,2,2,0.01,'AAT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,5 ,GAA.beta,GAA.zeta,2,2,0.01,'GAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,6 ,GAG.beta,GAG.zeta,2,2,0.01,'GAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,7 ,GAC.beta,GAC.zeta,2,2,0.01,'GAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,8 ,GAT.beta,GAT.zeta,2,2,0.01,'GAT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,9 ,CAA.beta,CAA.zeta,2,2,0.01,'CAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,10,CAG.beta,CAG.zeta,2,2,0.01,'CAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,11,CAC.beta,CAC.zeta,2,2,0.01,'CAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,12,CAT.beta,CAT.zeta,2,2,0.01,'CAT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,13,TAA.beta,TAA.zeta,2,2,0.01,'TAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,14,TAG.beta,TAG.zeta,2,2,0.01,'TAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,15,TAC.beta,TAC.zeta,2,2,0.01,'TAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,16,TAT.beta,TAT.zeta,2,2,0.01,'TAT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,17,AGA.beta,AGA.zeta,2,2,0.01,'AGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,18,AGG.beta,AGG.zeta,2,2,0.01,'AGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,19,AGC.beta,AGC.zeta,2,2,0.01,'AGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,20,AGT.beta,AGT.zeta,2,2,0.01,'AGT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,21,GGA.beta,GGA.zeta,2,2,0.01,'GGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,22,GGG.beta,GGG.zeta,2,2,0.01,'GGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,23,GGC.beta,GGC.zeta,2,2,0.01,'GGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,24,GGT.beta,GGT.zeta,2,2,0.01,'GGT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,25,CGA.beta,CGA.zeta,2,2,0.01,'CGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,26,CGG.beta,CGG.zeta,2,2,0.01,'CGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,27,CGC.beta,CGC.zeta,2,2,0.01,'CGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,28,CGT.beta,CGT.zeta,2,2,0.01,'CGT',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,29,TGA.beta,TGA.zeta,2,2,0.01,'TGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,30,TGG.beta,TGG.zeta,2,2,0.01,'TGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,31,TGC.beta,TGC.zeta,2,2,0.01,'TGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,32,TGT.beta,TGT.zeta,2,2,0.01,'TGT',' ',' ',0,60,0,60,colorbar=False)

pdf = matplotlib.backends.backend_pdf.PdfPages(name)
for fig in xrange(1,7):
    pdf.savefig(fig)

pdf.close()

pl.show()











