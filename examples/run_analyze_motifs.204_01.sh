#!/usr/bin/env bash

analyze_motifs \
--clusters inputs/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.Fc3Pv3.bed \
--species hg19_v19 \
--k 6 \
--nrand 1 \
--out_file outputs/204_01_RBFOX2.svg \
--pickle outputs/204_01_RBFOX2.pickle \
--genome_fasta inputs/all.fa \
--homer_out outputs/204_01_RBFOX2_homer_out
