#!/usr/bin/env bash

conda create -y -n analyze_motifs python=2.7

source activate analyze_motifs;

conda install -c anaconda -y \
matplotlib \
numpy \
pandas \
htseq \
seaborn \
scikit-learn;

conda install -c bioconda -y \
bedtools=2.26 pybedtools \
bx-python \
emboss \
homer;

# bedtools 2.27.1 breaks; BED6 files cause segfaults.
# unless homer/emboss is already in the path

pip install cairosvg==1.0.22; # cairo 2.0 is python3+ compatible only

## These are for C-libraries that compseq needs ##
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64/:$LD_LIBRARY_PATH
