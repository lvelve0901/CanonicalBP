#!/usr/bin/python

import pandas as pd
from ensembtool.base import *

abg = pd.read_table('Output_abg.txt',delim_whitespace=True)
seq = pd.read_table('Output_seq.txt',delim_whitespace=True)
abg_new = abg.loc[(abg.RMSD1 <= 2.0) & (abg.RMSD2 <= 2.0) & (seq.sf == 'B')]

print(len(abg_new))

gamma = abg_new.gamma.values*np.pi/180
beta = abg_new.beta.values
zeta = abg_new.zeta.values

fig = pl.figure(1,figsize=(12,8))

hist2D(fig,2,2,2,beta,zeta,1,1,0.004,' ',' ',' ',0,60,0,60,colorbar=True)

fig.subplots_adjust(right=0.8)

bins = [np.arange(-np.pi,np.pi+np.pi/30,np.pi/30),np.arange(0,60+2,2)]
H, x, y = np.histogram2d(gamma,beta,bins=bins,normed=1)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)    

ax = fig.add_subplot(221,projection='polar')
ax.pcolor(x,y,Hmasked)
ax.set_rticks([0, 15, 30, 45, 60])
ax.set_rlabel_position(-45)
ax.set_xticks(np.pi/180.*np.linspace(180,-180,8,endpoint=False))
ax.set_xticklabels([180,135,90,45,0,-45,-90,-135])
ax.grid(True)
ax.pcolor(x,y,Hmasked,norm=colors.Normalize(vmin=0,vmax=0.06))

bins = [np.arange(-np.pi,np.pi+np.pi/30,np.pi/30),np.arange(0,60+2,2)]
H, x, y = np.histogram2d(gamma,zeta,bins=bins,normed=1)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)    

ax = fig.add_subplot(224,projection='polar')
ax.pcolor(x,y,Hmasked)
ax.set_rticks([0, 15, 30, 45, 60])
ax.set_rlabel_position(-45)
ax.set_xticks(np.pi/180.*np.linspace(180,-180,8,endpoint=False))
ax.set_xticklabels([180,135,90,45,0,-45,-90,-135])
ax.grid(True)
ax.pcolor(x,y,Hmasked,norm=colors.Normalize(vmin=0,vmax=0.06))

pl.show()
