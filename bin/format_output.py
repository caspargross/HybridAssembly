#! /usr/bin/env python

import sys
import os
from Bio import SeqIO


# Assign variables from command line parameters

inputFile = sys.argv[1]
sampleID = sys.argv[2]
assemblyType = sys.argv[3]
minLength = int(sys.argv[4])
data_type = sys.argv[5]

long_contigs = []
input_handle=open(inputFile, 'rU')
output_handle=open(sampleID+'_'+assemblyType+'_genome.fasta', 'w')

for index, record in enumerate(SeqIO.parse(input_handle, 'fasta')):
    if len(record.seq) >= minLength:
        record.id = sampleID +'_'+ str(index+1)
        record.description = "assembler="+assemblyType+" length=" + str(len(record.seq)) + " datasource=" + data_type
        long_contigs.append(record)

SeqIO.write(long_contigs, output_handle, "fasta")

input_handle.close()
output_handle.close()
2

