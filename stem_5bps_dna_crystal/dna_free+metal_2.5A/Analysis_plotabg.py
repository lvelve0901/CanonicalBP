import pandas as pd
from ensembtool.base import *
import matplotlib.backends.backend_pdf

name = "Plotabg_dna_free+metal.pdf"
abg = pd.read_table('../Output_abg.txt',delim_whitespace=True)
seq = pd.read_table('../Output_seq.txt',delim_whitespace=True)
l = [i for s in read('dna_free+metal_2.5A.txt') for i in s]

abg_new = abg.loc[(abg.RMSD1 <= 2.0) & (abg.RMSD2 <= 2.0) & (seq.sf == 'B') & (seq.pdbid.isin(l))]


a = ABG(abg_new)
#a = ABG('Output_abg.txt')

print(len(abg_new))
print(len(abg_new.pdbid.unique()))
print "beta: " + str(np.mean(abg_new.beta)) + " " + str(np.std(abg_new.beta))
a.phl2()

pdf = matplotlib.backends.backend_pdf.PdfPages(name)
for fig in xrange(1,4):
    pdf.savefig(fig)

pdf.close()

pl.show()
