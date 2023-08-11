#!/usr/bin/env bash

#SBATCH --account=bgmp
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G

conda activate bgmp_py311

file_R1='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
file_R2='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
file_R3='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
file_R4='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'


/usr/bin/time -v ./demultiplexing.py -f1 $file_R1 -f2 $file_R2 -f3 $file_R3 -f4 $file_R4 -l 8 -a analysis9.txt -p matched_pairs.png