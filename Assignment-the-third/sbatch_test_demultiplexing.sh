#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1
#SBATCH --job-name=demult_test           #optional: job name
#SBATCH --output=demult_test_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demult_test_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

/usr/bin/time -v ./demultiplexing.py -r1 "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R1_input_test.fq" \
 -r2 "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R4_input_test.fq" \
 -I1 "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R2_input_test.fq" \
 -I2 "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R3_input_test.fq" \
 -EI "/projects/bgmp/shared/2017_sequencing/indexes.txt" -qc 29