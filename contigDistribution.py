# Small script to plot sequence length dsitrinution in a fasta file
# Caspar Gross 2018

import sys
from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

lengths = map(len, SeqIO.parse(sys.argv[1], 'fasta'))
print(lengths)

fig = pyplot.figure()
axes = pd.Series(lengths).hist(color='blue', bins = 1000)
fig = pyplot.gcf()
title = '' + sys.argv[2]
axes.set_title(title)
axes.set_xscale('symlog')
axes.set_xlabel("Nucleotides")
axes.set_ylabel("Number of Contigs")
axes.xaxis.grid(False)


ax = pd.Series(lengths).hist(color='blue', bins = 1000)
fig = ax.get_figure()
fig.savefig("dist" + sys.argv[2] + ".pdf", format='pdf')
