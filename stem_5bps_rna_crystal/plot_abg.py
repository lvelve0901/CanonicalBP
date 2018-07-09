import pandas as pd
from ensembtool.base import *
import matplotlib.backends.backend_pdf

abg = pd.read_table('Output_abg.txt',delim_whitespace=True)
seq = pd.read_table('Output_seq.txt',delim_whitespace=True)
abg_new = abg.loc[(abg.RMSD1 <= 2.0) & (abg.RMSD2 <= 2.0) & (seq.sf == 'A')]


a = ABG(abg_new)
#a = ABG('Output_abg.txt')

print(len(abg_new))
a.phl2()

pdf = matplotlib.backends.backend_pdf.PdfPages("rna_5bps.pdf")
for fig in xrange(1,4):
    pdf.savefig(fig)

pdf.close()

pl.show()

