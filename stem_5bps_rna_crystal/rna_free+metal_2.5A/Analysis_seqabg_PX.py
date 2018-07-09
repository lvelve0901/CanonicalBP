
#!/usr/bin/python

import pandas as pd
import pylab as pl
import numpy as np
import matplotlib.backends.backend_pdf
from commontool import read, plot2D, hist1D, hist2D, phist2D

name = "Plotseq_rna_free+metal_PX.pdf"

abg_init = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq_init = pd.read_table('../Output_seq.txt',delim_whitespace=True)
l = [i for s in read('rna_free+metal_2.5A.txt') for i in s]

abg = abg_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]
seq = seq_init.loc[(abg_init.RMSD1 <= 2.0) & (abg_init.RMSD2 <= 2.0) & (seq_init.sf == 'A') & (seq_init.pdbid.isin(l))]

AA = abg.loc[seq['3bps'].str[1:3] == 'AA']
AG = abg.loc[seq['3bps'].str[1:3] == 'AG']
AC = abg.loc[seq['3bps'].str[1:3] == 'AC']
AU = abg.loc[seq['3bps'].str[1:3] == 'AU']
GA = abg.loc[seq['3bps'].str[1:3] == 'GA']
GG = abg.loc[seq['3bps'].str[1:3] == 'GG']
GC = abg.loc[seq['3bps'].str[1:3] == 'GC']
GU = abg.loc[seq['3bps'].str[1:3] == 'GU']

fig1 = pl.figure(1,figsize=(8,13))

fig1.suptitle('beta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.beta,1,'AA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,2, AG.beta,1,'AG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,3, AC.beta,1,'AC','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,4, AU.beta,1,'AU','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,5, GA.beta,1,'GA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,6, GG.beta,1,'GG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,7, GC.beta,1,'GC','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,8, GU.beta,1,'GU','red',0,60,0.2,legend=None)

fig1 = pl.figure(2,figsize=(8,13))

fig1.suptitle('gamma distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.gamma,6,'AA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,2, AG.gamma,6,'AG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,3, AC.gamma,6,'AC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,4, AU.gamma,6,'AU','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,5, GA.gamma,6,'GA','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,6, GG.gamma,6,'GG','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,7, GC.gamma,6,'GC','red',-180,180,0.1,legend=None)
hist1D(fig1,4,2,8, GU.gamma,6,'GU','red',-180,180,0.1,legend=None)

fig1 = pl.figure(3,figsize=(8,13))

fig1.suptitle('zeta distribution',fontsize=14,fontweight='bold')

hist1D(fig1,4,2,1, AA.zeta,1,'AA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,2, AG.zeta,1,'AG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,3, AC.zeta,1,'AC','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,4, AU.zeta,1,'AU','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,5, GA.zeta,1,'GA','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,6, GG.zeta,1,'GG','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,7, GC.zeta,1,'GC','red',0,60,0.2,legend=None)
hist1D(fig1,4,2,8, GU.zeta,1,'GU','red',0,60,0.2,legend=None)

fig1 = pl.figure(4,figsize=(8,13))

fig1.suptitle('gamma-beta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,2,1,AA.gamma*np.pi/180,AA.beta,np.pi/20,2,0.04,'AA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,2,AG.gamma*np.pi/180,AG.beta,np.pi/20,2,0.04,'AG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,3,AC.gamma*np.pi/180,AC.beta,np.pi/20,2,0.04,'AC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,4,AU.gamma*np.pi/180,AU.beta,np.pi/20,2,0.04,'AU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,5,GA.gamma*np.pi/180,GA.beta,np.pi/20,2,0.04,'GA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,6,GG.gamma*np.pi/180,GG.beta,np.pi/20,2,0.04,'GG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,7,GC.gamma*np.pi/180,GC.beta,np.pi/20,2,0.04,'GC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,8,GU.gamma*np.pi/180,GU.beta,np.pi/20,2,0.04,'GU',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(5,figsize=(8,13))

fig1.suptitle('gamma-zeta distribution',fontsize=14,fontweight='bold')

phist2D(fig1,4,2,1,AA.gamma*np.pi/180,AA.zeta,np.pi/20,2,0.04,'AA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,2,AG.gamma*np.pi/180,AG.zeta,np.pi/20,2,0.04,'AG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,3,AC.gamma*np.pi/180,AC.zeta,np.pi/20,2,0.04,'AC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,4,AU.gamma*np.pi/180,AU.zeta,np.pi/20,2,0.04,'AU',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,5,GA.gamma*np.pi/180,GA.zeta,np.pi/20,2,0.04,'GA',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,6,GG.gamma*np.pi/180,GG.zeta,np.pi/20,2,0.04,'GG',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,7,GC.gamma*np.pi/180,GC.zeta,np.pi/20,2,0.04,'GC',0,60,[0,15,30,45,60],'')
phist2D(fig1,4,2,8,GU.gamma*np.pi/180,GU.zeta,np.pi/20,2,0.04,'GU',0,60,[0,15,30,45,60],'')

fig1 = pl.figure(6,figsize=(8,13))

fig1.suptitle('beta-zeta distribution',fontsize=14,fontweight='bold')

hist2D(fig1,4,2,1,AA.beta,AA.zeta,2,2,0.006,'AA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,2,AG.beta,AG.zeta,2,2,0.006,'AG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,3,AC.beta,AC.zeta,2,2,0.006,'AC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,4,AU.beta,AU.zeta,2,2,0.006,'AU',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,5,GA.beta,GA.zeta,2,2,0.006,'GA',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,6,GG.beta,GG.zeta,2,2,0.006,'GG',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,7,GC.beta,GC.zeta,2,2,0.006,'GC',' ',' ',0,60,0,60,colorbar=False)
hist2D(fig1,4,2,8,GU.beta,GU.zeta,2,2,0.006,'GU',' ',' ',0,60,0,60,colorbar=False)

pdf = matplotlib.backends.backend_pdf.PdfPages(name)
for fig in xrange(1,7):
    pdf.savefig(fig)

pdf.close()

pl.show()


