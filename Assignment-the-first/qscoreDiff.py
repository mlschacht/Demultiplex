#!/usr/bin/env python

import argparse
import bioinfo
import matplotlib.pyplot as plt
import gzip


#Un-comment when ready for the real data
Read1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
Read2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
Index1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
Index2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"

#Comment out when ready for real data
# Read1 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R1_input_test.fq"
# Read2 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R4_input_test.fq"
# Index1 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R2_input_test.fq"
# Index2 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R3_input_test.fq"

Read1_qscores, Read1_lines = bioinfo.populate_list_gzip(Read1)
Read2_qscores, Read2_lines = bioinfo.populate_list_gzip(Read2)
Index1_qscores, Index1_lines = bioinfo.populate_list_gzip(Index1) #list is still 101 long...
Index2_qscores, Index2_lines = bioinfo.populate_list_gzip(Index2) #list is still 101 long...

#update the qscores list to be the average at each position
for pos, q in enumerate(Read1_qscores):
    Read1_qscores[pos] = Read1_qscores[pos] / (Read1_lines/4)
    Read2_qscores[pos] = Read2_qscores[pos] / (Read2_lines/4)
    Index1_qscores[pos] = Index1_qscores[pos] / (Index1_lines/4)
    Index2_qscores[pos] = Index2_qscores[pos] / (Index2_lines/4)

#Read 1 distribution Graph
plt.bar(range(0, 101), Read1_qscores, width=0.8)
plt.title('Average Read 1 Quality Score for Each Base Pair')
plt.xlabel('Base Pair Position on Sequence')
plt.ylabel('Average Phred 33 Quality Score')
plt.savefig("Read1_dist.png")
plt.cla()

#Read 2 distribution Graph
plt.bar(range(0, 101), Read2_qscores, width=0.8)
plt.title('Average Read 2 Quality Score for Each Base Pair')
plt.xlabel('Base Pair Position on Sequence')
plt.ylabel('Average Phred 33 Quality Score')
plt.savefig("Read2_dist.png")
plt.cla()

#Index 1 distribution Graph
plt.bar(range(0, 8), Index1_qscores[:8], width=0.8)
plt.title('Average Index 1 Quality Score for Each Base Pair')
plt.xlabel('Base Pair Position on Sequence')
plt.ylabel('Average Phred 33 Quality Score')
plt.savefig("Index1_dist.png")
plt.cla()

#Index 2 distribution Graph
plt.bar(range(0, 8), Index2_qscores[:8], width=0.8)
plt.title('Average Index 2 Quality Score for Each Base Pair')
plt.xlabel('Base Pair Position on Sequence')
plt.ylabel('Average Phred 33 Quality Score')
plt.savefig("Index2_dist.png")
plt.cla()