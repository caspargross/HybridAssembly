#!/usr/bin/env python

# Calculate summary stats and create plots

import sys
import os
import json
import codecs
from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle, islice
import scipy
import re
import xml.etree.ElementTree as ET
import numpy as np


# Get environment variables ---------------------------------------------------
id = sys.argv[1].strip()
assemblers = sys.argv[2].strip('][').strip().split(',')
assemblers = [x.strip(' ') for x in assemblers]
genomes = sys.argv[3].strip('][').strip().split(',')
genomes = [x.strip(' ') for x in genomes]
readstats = sys.argv[4].strip('][').strip().split(',')
readstats = [x.strip(' ') for x in readstats]
stattypes = sys.argv[5].strip('][').strip().split(',')
stattypes = [x.strip(' ') for x in stattypes]

print(id)
print(assemblers)
print(genomes)
print(readstats)
print(stattypes)
print(genomes)



# Parse NanoPlot ---------------------------------------------------------------
def parseNanoPlot(npFile):

    SECTIONS="summary|cutoff|highest mean|longest read"
    dat_np = {}

    with open(npFile) as f:
        lines = f.readlines()
        grp = ""

    for line in lines:

        m = re.search(SECTIONS, line)
        if m is not None:
            grp = m.group()
            dat_np[grp] = {}

        else:
            line = line.strip()
            #print("line: " + line)
            # Split on tabs and two or more white spaces
            lsp = re.split(r'\t|\s{2,}', line)
            data_id = re.sub('[^A-Za-z0-9]+', '', lsp[0].lower())
            dat_np[grp][data_id] = {}
            dat_np[grp][data_id]["data_description"] = lsp[0]
            dat_np[grp][data_id]["value_num"] = float(lsp[1].split(' ')[0].replace(',',""))
            dat_np[grp][data_id]["value"] = lsp[1]

    return(dat_np)

# Parse seqPurge ---------------------------------------------------------------
def parseSeqPurge(spFile):

    xroot = ET.parse(spFile).getroot()[0]
    xqualPar = xroot.findall('{http://www.prime-xs.eu/ms/qcml}qualityParameter')
    d_seqpurge = {}

    for child in xqualPar:
        data_id = re.sub('[^A-Za-z0-9]+', '', child.attrib["name"].lower())
        d_seqpurge[data_id] = {}
        d_seqpurge[data_id]["data_description"] = child.attrib["description"]
        d_seqpurge[data_id]["value_num"] = float(re.search('[0-9.,]*$',child.attrib["value"]).group(0))
        d_seqpurge[data_id]["value"] = child.attrib["value"]

    return(d_seqpurge)

# Calculate Assembly stats ---------------------------------------------------
def calcAssemblyStats(assemblyFasta):

    d_assemb = {}
    seqLengths = []
    totalLength = 0
    nContigs = 0

    for record in (SeqIO.parse(assemblyFasta, 'fasta')):
        bp = len(record.seq)
        seqLengths.append(bp)
        totalLength += bp
        nContigs += 1

    seqLengths = np.array(sorted(seqLengths))
    seqLengths = seqLengths[::-1]
    seqLengthsSums = np.cumsum(seqLengths).astype(int)

    N50 = int(min(seqLengthsSums[seqLengthsSums > (totalLength * 0.5)]))
    N90 = int(min(seqLengthsSums[seqLengthsSums > (totalLength * 0.9)]))

    L50 = int(np.where(seqLengthsSums==N50)[0])
    L90 = int(np.where(seqLengthsSums == N90)[0])


    d_assemb["N50"] = N50
    d_assemb["N90"] = N90
    d_assemb["L50"] = L50
    d_assemb["L90"] = L90
    d_assemb["n_contigs"] = nContigs
    d_assemb["total_length"] = totalLength
    d_assemb["lengths"] = seqLengths.tolist()

    nx_x = (seqLengthsSums / totalLength)
    nx_x = np.insert(nx_x, 0, 0)
    nx_x = np.append(nx_x, 1)
    nx_x.astype(float)
    nx_y = seqLengths
    nx_y = np.insert(nx_y, 0, max(seqLengths))
    nx_y = np.append(nx_y, 0)
    nx_y.astype(int)

    d_assemb["nx_x"] = nx_x.tolist()
    d_assemb["nx_y"] = nx_y.tolist()

    return(d_assemb)

# Read and Import data -----------------------------------------------

dat = {id:{el:{} for el in assemblers}}

for i, a in enumerate(assemblers):
    dat[id][a]=calcAssemblyStats(genomes[i])

for i, s in enumerate(stattypes):
    if s == "read_qc":
        dat[id][s] = parseSeqPurge(readstats[i])
    else:
        dat[id][s] = parseNanoPlot(readstats[i])

# Calculate coverages
for a in assemblers:
    dat[id][a]["lr_cov_raw"] = dat[id]["raw"]["summary"]["totalbases"]["value_num"] / dat[id][a]["total_length"]
    dat[id][a]["lr_cov_filtered"] = dat[id]["filtered"]["summary"]["totalbases"]["value_num"] / dat[id][a]["total_length"]
    if "read_qc" in dat[id]:
        dat[id][a]["sr_cov"] = (dat[id]["read_qc"]["basessequencedmb"]["value_num"] * 1e6)/ dat[id][a]["total_length"]

# Create plots -------------------------------------------------------

# NX_stepplot
for a in assemblers:
    plt.step(dat[id][a]["nx_x"], dat[id][a]["nx_y"], where='post')
plt.xlabel("Nx")
plt.ylabel("Contig length (kbp)")
plt.legend(assemblers)
plt.title("NX plot")
plt.savefig(fname="nx_plot.pdf")
plt.savefig(fname="nx_plot.png")

# Contig Lengths plot
lengths = [dat[id][x]["lengths"] for x in assemblers]
df_len = pd.DataFrame([dat[id][x]["lengths"] for x in assemblers], index=assemblers)
N50s = [dat[id][x]["N50"] for x in assemblers]
N90s = [dat[id][x]["N50"] for x in assemblers]

#lengths.rows = assemblers
print(lengths)
print(df_len)
print(N50s)
barColors = list(islice(cycle(['chartreuse', 'skyblue']), None, len(df_len)))
fig, ax = plt.subplots()
df_len.plot.barh(stacked=True,
                      legend=False,
                      color = barColors,
                      edgecolor = "k", ax = ax)
ax.set_xlabel("Cumulative contig length")
ax.set_ylabel("Assembly types")
plt.savefig(fname="contigs.pdf")
plt.savefig(fname="contigs.png")
#plt.show()


# Export data:
with open('qc_summary_' + id + '.json', 'w') as fp:
    json.dump(dat, fp)
