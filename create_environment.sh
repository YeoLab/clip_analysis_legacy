#!/usr/bin/env bash

conda create -n analyze_motifs python=2.7

source activate analyze_motifs;

conda install -y matplotlib numpy pandas htseq seaborn;
conda install -y bedtools=2.26 pybedtools; # bedtools 2.27.1 breaks; BED6 files cause segfaults.
pip install cairosvg==1.0.22; # cairo 2.0 is python3+ compatible only
conda install -y bx-python;
conda install -y scikit-learn;
conda install -y -c bioconda emboss; # unless emboss is already in the path
conda install -y -c bioconda homer; # unless homer is already in the path
# already in package directory
python setup.py install