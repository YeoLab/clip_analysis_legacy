#!/usr/bin/env bash

python analyze_motifs.py \
--clusters /home/bay001/projects/parp13_ago2_20171015/permanent_data/eCLIP-0.1.7/05_clip_analysis/inputs/PARP13.KO3P1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed \
--species hg19 \
--k 6 \
--nrand 1 \
--out_file /projects/ps-yeolab3/bay001/tmp/out_file.svg \
--pickle /projects/ps-yeolab3/bay001/tmp/pickle.pickle \
--genome_fasta /projects/ps-yeolab/genomes/hg19/chromosomes/all.fa \
--homer_out /projects/ps-yeolab3/bay001/tmp/homer_out/