
#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
import matplotlib.backends.backend_pdf
from commontool import read, plot2D, hist1D, hist2D, phist2D

name = "Plotseq_rna_free+metal_XPX.pdf"

abg_init = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('../Output_seq.txt',delim_whitespace=True)
l = [i for s in read('rna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]

AAA = abg.loc[seq['3bps'] == 'AAA']
AAG = abg.loc[seq['3bps'] == 'AAG']
AAC = abg.loc[seq['3bps'] == 'AAC']
AAU = abg.loc[seq['3bps'] == 'AAU']
GAA = abg.loc[seq['3bps'] == 'GAA']
GAG = abg.loc[seq['3bps'] == 'GAG']
GAC = abg.loc[seq['3bps'] == 'GAC']
GAU = abg.loc[seq['3bps'] == 'GAU']
CAA = abg.loc[seq['3bps'] == 'CAA']
CAG = abg.loc[seq['3bps'] == 'CAG']
CAC = abg.loc[seq['3bps'] == 'CAC']
CAU = abg.loc[seq['3bps'] == 'CAU']
UAA = abg.loc[seq['3bps'] == 'UAA']
UAG = abg.loc[seq['3bps'] == 'UAG']
UAC = abg.loc[seq['3bps'] == 'UAC']
UAU = abg.loc[seq['3bps'] == 'UAU']
AGA = abg.loc[seq['3bps'] == 'AGA']
AGG = abg.loc[seq['3bps'] == 'AGG']
AGC = abg.loc[seq['3bps'] == 'AGC']
AGU = abg.loc[seq['3bps'] == 'AGU']
GGA = abg.loc[seq['3bps'] == 'GGA']
GGG = abg.loc[seq['3bps'] == 'GGG']
GGC = abg.loc[seq['3bps'] == 'GGC']
GGU = abg.loc[seq['3bps'] == 'GGU']
CGA = abg.loc[seq['3bps'] == 'CGA']
CGG = abg.loc[seq['3bps'] == 'CGG']
CGC = abg.loc[seq['3bps'] == 'CGC']
CGU = abg.loc[seq['3bps'] == 'CGU']
UGA = abg.loc[seq['3bps'] == 'UGA']
UGG = abg.loc[seq['3bps'] == 'UGG']
UGC = abg.loc[seq['3bps'] == 'UGC']
UGU = abg.loc[seq['3bps'] == 'UGU']

fig1 = pl.figure(1,figsize=(27,13))

fig1.suptitle('beta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.beta,1,'AAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,2 , AAG.beta,1,'AAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,3 , AAC.beta,1,'AAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,4 , AAU.beta,1,'AAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,5 , GAA.beta,1,'GAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,6 , GAG.beta,1,'GAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,7 , GAC.beta,1,'GAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,8 , GAU.beta,1,'GAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,9 , CAA.beta,1,'CAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,10, CAG.beta,1,'CAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,11, CAC.beta,1,'CAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,12, CAU.beta,1,'CAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,13, UAA.beta,1,'UAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,14, UAG.beta,1,'UAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,15, UAC.beta,1,'UAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,16, UAU.beta,1,'UAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,17, AGA.beta,1,'AGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,18, AGG.beta,1,'AGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,19, AGC.beta,1,'AGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,20, AGU.beta,1,'AGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,21, GGA.beta,1,'GGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,22, GGG.beta,1,'GGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,23, GGC.beta,1,'GGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,24, GGU.beta,1,'GGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,25, CGA.beta,1,'CGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,26, CGG.beta,1,'CGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,27, CGC.beta,1,'CGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,28, CGU.beta,1,'CGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,29, UGA.beta,1,'UGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,30, UGG.beta,1,'UGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,31, UGC.beta,1,'UGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,32, UGU.beta,1,'UGU','red',0,60,0.2,legend=None)

fig1 = pl.figure(2,figsize=(27,13))

fig1.suptitle('gamma distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.gamma,6,'AAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,2 , AAG.gamma,6,'AAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,3 , AAC.gamma,6,'AAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,4 , AAU.gamma,6,'AAU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,5 , GAA.gamma,6,'GAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,6 , GAG.gamma,6,'GAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,7 , GAC.gamma,6,'GAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,8 , GAU.gamma,6,'GAU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,9 , CAA.gamma,6,'CAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,10, CAG.gamma,6,'CAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,11, CAC.gamma,6,'CAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,12, CAU.gamma,6,'CAU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,13, UAA.gamma,6,'UAA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,14, UAG.gamma,6,'UAG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,15, UAC.gamma,6,'UAC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,16, UAU.gamma,6,'UAU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,17, AGA.gamma,6,'AGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,18, AGG.gamma,6,'AGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,19, AGC.gamma,6,'AGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,20, AGU.gamma,6,'AGU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,21, GGA.gamma,6,'GGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,22, GGG.gamma,6,'GGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,23, GGC.gamma,6,'GGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,24, GGU.gamma,6,'GGU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,25, CGA.gamma,6,'CGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,26, CGG.gamma,6,'CGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,27, CGC.gamma,6,'CGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,28, CGU.gamma,6,'CGU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,29, UGA.gamma,6,'UGA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,30, UGG.gamma,6,'UGG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,31, UGC.gamma,6,'UGC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,8,32, UGU.gamma,6,'UGU','red',-180,180,0.1,legend=None)

fig1 = pl.figure(3,figsize=(27,13))

fig1.suptitle('zeta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,8,1 , AAA.zeta,1,'AAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,2 , AAG.zeta,1,'AAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,3 , AAC.zeta,1,'AAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,4 , AAU.zeta,1,'AAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,5 , GAA.zeta,1,'GAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,6 , GAG.zeta,1,'GAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,7 , GAC.zeta,1,'GAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,8 , GAU.zeta,1,'GAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,9 , CAA.zeta,1,'CAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,10, CAG.zeta,1,'CAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,11, CAC.zeta,1,'CAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,12, CAU.zeta,1,'CAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,13, UAA.zeta,1,'UAA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,14, UAG.zeta,1,'UAG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,15, UAC.zeta,1,'UAC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,16, UAU.zeta,1,'UAU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,17, AGA.zeta,1,'AGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,18, AGG.zeta,1,'AGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,19, AGC.zeta,1,'AGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,20, AGU.zeta,1,'AGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,21, GGA.zeta,1,'GGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,22, GGG.zeta,1,'GGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,23, GGC.zeta,1,'GGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,24, GGU.zeta,1,'GGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,25, CGA.zeta,1,'CGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,26, CGG.zeta,1,'CGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,27, CGC.zeta,1,'CGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,28, CGU.zeta,1,'CGU','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,29, UGA.zeta,1,'UGA','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,30, UGG.zeta,1,'UGG','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,31, UGC.zeta,1,'UGC','red',0,60,0.2,legend=None)
hist1D(fig1,4,8,32, UGU.zeta,1,'UGU','red',0,60,0.2,legend=None)

fig1 = pl.figure(4,figsize=(27,13))

fig1.suptitle('gamma-beta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,8,1 ,AAA.gamma*np.pi/180,AAA.beta,np.pi/20,2,0.08,'AAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,2 ,AAG.gamma*np.pi/180,AAG.beta,np.pi/20,2,0.08,'AAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,3 ,AAC.gamma*np.pi/180,AAC.beta,np.pi/20,2,0.08,'AAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,4 ,AAU.gamma*np.pi/180,AAU.beta,np.pi/20,2,0.08,'AAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,5 ,GAA.gamma*np.pi/180,GAA.beta,np.pi/20,2,0.08,'GAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,6 ,GAG.gamma*np.pi/180,GAG.beta,np.pi/20,2,0.08,'GAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,7 ,GAC.gamma*np.pi/180,GAC.beta,np.pi/20,2,0.08,'GAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,8 ,GAU.gamma*np.pi/180,GAU.beta,np.pi/20,2,0.08,'GAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,9 ,CAA.gamma*np.pi/180,CAA.beta,np.pi/20,2,0.08,'CAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,10,CAG.gamma*np.pi/180,CAG.beta,np.pi/20,2,0.08,'CAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,11,CAC.gamma*np.pi/180,CAC.beta,np.pi/20,2,0.08,'CAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,12,CAU.gamma*np.pi/180,CAU.beta,np.pi/20,2,0.08,'CAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,13,UAA.gamma*np.pi/180,UAA.beta,np.pi/20,2,0.08,'UAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,14,UAG.gamma*np.pi/180,UAG.beta,np.pi/20,2,0.08,'UAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,15,UAC.gamma*np.pi/180,UAC.beta,np.pi/20,2,0.08,'UAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,16,UAU.gamma*np.pi/180,UAU.beta,np.pi/20,2,0.08,'UAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,17,AGA.gamma*np.pi/180,AGA.beta,np.pi/20,2,0.08,'AGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,18,AGG.gamma*np.pi/180,AGG.beta,np.pi/20,2,0.08,'AGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,19,AGC.gamma*np.pi/180,AGC.beta,np.pi/20,2,0.08,'AGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,20,AGU.gamma*np.pi/180,AGU.beta,np.pi/20,2,0.08,'AGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,21,GGA.gamma*np.pi/180,GGA.beta,np.pi/20,2,0.08,'GGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,22,GGG.gamma*np.pi/180,GGG.beta,np.pi/20,2,0.08,'GGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,23,GGC.gamma*np.pi/180,GGC.beta,np.pi/20,2,0.08,'GGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,24,GGU.gamma*np.pi/180,GGU.beta,np.pi/20,2,0.08,'GGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,25,CGA.gamma*np.pi/180,CGA.beta,np.pi/20,2,0.08,'CGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,26,CGG.gamma*np.pi/180,CGG.beta,np.pi/20,2,0.08,'CGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,27,CGC.gamma*np.pi/180,CGC.beta,np.pi/20,2,0.08,'CGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,28,CGU.gamma*np.pi/180,CGU.beta,np.pi/20,2,0.08,'CGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,29,UGA.gamma*np.pi/180,UGA.beta,np.pi/20,2,0.08,'UGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,30,UGG.gamma*np.pi/180,UGG.beta,np.pi/20,2,0.08,'UGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,31,UGC.gamma*np.pi/180,UGC.beta,np.pi/20,2,0.08,'UGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,32,UGU.gamma*np.pi/180,UGU.beta,np.pi/20,2,0.08,'UGU',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(5,figsize=(27,13))

fig1.suptitle('gamma-zeta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,8,1 ,AAA.gamma*np.pi/180,AAA.zeta,np.pi/20,2,0.08,'AAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,2 ,AAG.gamma*np.pi/180,AAG.zeta,np.pi/20,2,0.08,'AAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,3 ,AAC.gamma*np.pi/180,AAC.zeta,np.pi/20,2,0.08,'AAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,4 ,AAU.gamma*np.pi/180,AAU.zeta,np.pi/20,2,0.08,'AAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,5 ,GAA.gamma*np.pi/180,GAA.zeta,np.pi/20,2,0.08,'GAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,6 ,GAG.gamma*np.pi/180,GAG.zeta,np.pi/20,2,0.08,'GAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,7 ,GAC.gamma*np.pi/180,GAC.zeta,np.pi/20,2,0.08,'GAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,8 ,GAU.gamma*np.pi/180,GAU.zeta,np.pi/20,2,0.08,'GAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,9 ,CAA.gamma*np.pi/180,CAA.zeta,np.pi/20,2,0.08,'CAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,10,CAG.gamma*np.pi/180,CAG.zeta,np.pi/20,2,0.08,'CAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,11,CAC.gamma*np.pi/180,CAC.zeta,np.pi/20,2,0.08,'CAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,12,CAU.gamma*np.pi/180,CAU.zeta,np.pi/20,2,0.08,'CAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,13,UAA.gamma*np.pi/180,UAA.zeta,np.pi/20,2,0.08,'UAA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,14,UAG.gamma*np.pi/180,UAG.zeta,np.pi/20,2,0.08,'UAG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,15,UAC.gamma*np.pi/180,UAC.zeta,np.pi/20,2,0.08,'UAC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,16,UAU.gamma*np.pi/180,UAU.zeta,np.pi/20,2,0.08,'UAU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,17,AGA.gamma*np.pi/180,AGA.zeta,np.pi/20,2,0.08,'AGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,18,AGG.gamma*np.pi/180,AGG.zeta,np.pi/20,2,0.08,'AGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,19,AGC.gamma*np.pi/180,AGC.zeta,np.pi/20,2,0.08,'AGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,20,AGU.gamma*np.pi/180,AGU.zeta,np.pi/20,2,0.08,'AGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,21,GGA.gamma*np.pi/180,GGA.zeta,np.pi/20,2,0.08,'GGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,22,GGG.gamma*np.pi/180,GGG.zeta,np.pi/20,2,0.08,'GGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,23,GGC.gamma*np.pi/180,GGC.zeta,np.pi/20,2,0.08,'GGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,24,GGU.gamma*np.pi/180,GGU.zeta,np.pi/20,2,0.08,'GGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,25,CGA.gamma*np.pi/180,CGA.zeta,np.pi/20,2,0.08,'CGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,26,CGG.gamma*np.pi/180,CGG.zeta,np.pi/20,2,0.08,'CGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,27,CGC.gamma*np.pi/180,CGC.zeta,np.pi/20,2,0.08,'CGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,28,CGU.gamma*np.pi/180,CGU.zeta,np.pi/20,2,0.08,'CGU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,29,UGA.gamma*np.pi/180,UGA.zeta,np.pi/20,2,0.08,'UGA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,30,UGG.gamma*np.pi/180,UGG.zeta,np.pi/20,2,0.08,'UGG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,31,UGC.gamma*np.pi/180,UGC.zeta,np.pi/20,2,0.08,'UGC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,8,32,UGU.gamma*np.pi/180,UGU.zeta,np.pi/20,2,0.08,'UGU',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(6,figsize=(27,13))

fig1.suptitle('beta-zeta distribution',fontsize=14,fontweight='bold')

hist2D(fig1,4,8,1 ,AAA.beta,AAA.zeta,2,2,0.01,'AAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,2 ,AAG.beta,AAG.zeta,2,2,0.01,'AAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,3 ,AAC.beta,AAC.zeta,2,2,0.01,'AAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,4 ,AAU.beta,AAU.zeta,2,2,0.01,'AAU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,5 ,GAA.beta,GAA.zeta,2,2,0.01,'GAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,6 ,GAG.beta,GAG.zeta,2,2,0.01,'GAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,7 ,GAC.beta,GAC.zeta,2,2,0.01,'GAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,8 ,GAU.beta,GAU.zeta,2,2,0.01,'GAU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,9 ,CAA.beta,CAA.zeta,2,2,0.01,'CAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,10,CAG.beta,CAG.zeta,2,2,0.01,'CAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,11,CAC.beta,CAC.zeta,2,2,0.01,'CAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,12,CAU.beta,CAU.zeta,2,2,0.01,'CAU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,13,UAA.beta,UAA.zeta,2,2,0.01,'UAA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,14,UAG.beta,UAG.zeta,2,2,0.01,'UAG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,15,UAC.beta,UAC.zeta,2,2,0.01,'UAC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,16,UAU.beta,UAU.zeta,2,2,0.01,'UAU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,17,AGA.beta,AGA.zeta,2,2,0.01,'AGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,18,AGG.beta,AGG.zeta,2,2,0.01,'AGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,19,AGC.beta,AGC.zeta,2,2,0.01,'AGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,20,AGU.beta,AGU.zeta,2,2,0.01,'AGU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,21,GGA.beta,GGA.zeta,2,2,0.01,'GGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,22,GGG.beta,GGG.zeta,2,2,0.01,'GGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,23,GGC.beta,GGC.zeta,2,2,0.01,'GGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,24,GGU.beta,GGU.zeta,2,2,0.01,'GGU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,25,CGA.beta,CGA.zeta,2,2,0.01,'CGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,26,CGG.beta,CGG.zeta,2,2,0.01,'CGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,27,CGC.beta,CGC.zeta,2,2,0.01,'CGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,28,CGU.beta,CGU.zeta,2,2,0.01,'CGU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,29,UGA.beta,UGA.zeta,2,2,0.01,'UGA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,30,UGG.beta,UGG.zeta,2,2,0.01,'UGG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,31,UGC.beta,UGC.zeta,2,2,0.01,'UGC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,8,32,UGU.beta,UGU.zeta,2,2,0.01,'UGU',' ',' ',0,60,0,60,colorbar=False)

pdf = matplotlib.backends.backend_pdf.PdfPages(name)
for fig in xrange(1,7):
    pdf.savefig(fig)

pdf.close()

pl.show()











