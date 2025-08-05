#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1
#SBATCH --job-name=demult           #optional: job name
#SBATCH --output=demult_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demult_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

/usr/bin/time -v ./demultiplexing.py -r1 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" \
 -r2 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"\
 -I1 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"\
 -I2 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"\
 -EI "/projects/bgmp/shared/2017_sequencing/indexes.txt"\
 -qc 29