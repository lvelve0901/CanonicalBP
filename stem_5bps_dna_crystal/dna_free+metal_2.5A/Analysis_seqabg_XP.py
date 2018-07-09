
#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
import matplotlib.backends.backend_pdf
from commontool import read, plot2D, hist1D, hist2D, phist2D

name = "Plotseq_dna_free+metal_XP.pdf"

abg_init = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('../Output_seq.txt',delim_whitespace=True)
l = [i for s in read('dna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'B') & (seq_init.pdbid.isin(l))]

AA = abg.loc[seq['3bps'].str[0:2] == 'AA']
AG = abg.loc[seq['3bps'].str[0:2] == 'AG']
GA = abg.loc[seq['3bps'].str[0:2] == 'GA']
GG = abg.loc[seq['3bps'].str[0:2] == 'GG']
CA = abg.loc[seq['3bps'].str[0:2] == 'CA']
CG = abg.loc[seq['3bps'].str[0:2] == 'CG']
TA = abg.loc[seq['3bps'].str[0:2] == 'TA']
TG = abg.loc[seq['3bps'].str[0:2] == 'TG']

fig1 = pl.figure(1,figsize=(8,13))

fig1.suptitle('beta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.beta,1,'AA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,2, AG.beta,1,'AG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,3, GA.beta,1,'GA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,4, GG.beta,1,'GG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,5, CA.beta,1,'CA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,6, CG.beta,1,'CG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,7, TA.beta,1,'TA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,8, TG.beta,1,'TG','red',0,60,0.2,legend=None)

fig1 = pl.figure(2,figsize=(8,13))

fig1.suptitle('gamma distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.gamma,6,'AA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,2, AG.gamma,6,'AG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,3, GA.gamma,6,'GA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,4, GG.gamma,6,'GG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,5, CA.gamma,6,'CA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,6, CG.gamma,6,'CG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,7, TA.gamma,6,'TA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,8, TG.gamma,6,'TG','red',-180,180,0.1,legend=None)

fig1 = pl.figure(3,figsize=(8,13))

fig1.suptitle('zeta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.zeta,1,'AA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,2, AG.zeta,1,'AG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,3, GA.zeta,1,'GA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,4, GG.zeta,1,'GG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,5, CA.zeta,1,'CA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,6, CG.zeta,1,'CG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,7, TA.zeta,1,'TA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,8, TG.zeta,1,'TG','red',0,60,0.2,legend=None)

fig1 = pl.figure(4,figsize=(8,13))

fig1.suptitle('gamma-beta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,2,1,AA.gamma*np.pi/180,AA.beta,np.pi/20,2,0.04,'AA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,2,AG.gamma*np.pi/180,AG.beta,np.pi/20,2,0.04,'AG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,3,GA.gamma*np.pi/180,GA.beta,np.pi/20,2,0.04,'GA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,4,GG.gamma*np.pi/180,GG.beta,np.pi/20,2,0.04,'GG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,5,CA.gamma*np.pi/180,CA.beta,np.pi/20,2,0.04,'CA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,6,CG.gamma*np.pi/180,CG.beta,np.pi/20,2,0.04,'CG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,7,TA.gamma*np.pi/180,TA.beta,np.pi/20,2,0.04,'TA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,8,TG.gamma*np.pi/180,TG.beta,np.pi/20,2,0.04,'TG',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(5,figsize=(8,13))

fig1.suptitle('gamma-zeta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,2,1,AA.gamma*np.pi/180,AA.zeta,np.pi/20,2,0.04,'AA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,2,AG.gamma*np.pi/180,AG.zeta,np.pi/20,2,0.04,'AG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,3,GA.gamma*np.pi/180,GA.zeta,np.pi/20,2,0.04,'GA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,4,GG.gamma*np.pi/180,GG.zeta,np.pi/20,2,0.04,'GG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,5,CA.gamma*np.pi/180,CA.zeta,np.pi/20,2,0.04,'CA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,6,CG.gamma*np.pi/180,CG.zeta,np.pi/20,2,0.04,'CG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,7,TA.gamma*np.pi/180,TA.zeta,np.pi/20,2,0.04,'TA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,8,TG.gamma*np.pi/180,TG.zeta,np.pi/20,2,0.04,'TG',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(6,figsize=(8,13))

fig1.suptitle('beta-zeta distribution',fontsize=14,fontweight='bold')

hist2D(fig1,4,2,1,AA.beta,AA.zeta,2,2,0.006,'AA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,2,AG.beta,AG.zeta,2,2,0.006,'AG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,3,GA.beta,GA.zeta,2,2,0.006,'GA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,4,GG.beta,GG.zeta,2,2,0.006,'GG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,5,CA.beta,CA.zeta,2,2,0.006,'CA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,6,CG.beta,CG.zeta,2,2,0.006,'CG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,7,TA.beta,TA.zeta,2,2,0.006,'TA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,8,TG.beta,TG.zeta,2,2,0.006,'TG',' ',' ',0,60,0,60,colorbar=False)

pdf = matplotlib.backends.backend_pdf.PdfPages(name)
for fig in xrange(1,7):
    pdf.savefig(fig)

pdf.close()

pl.show()











