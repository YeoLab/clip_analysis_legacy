# CLIP ANALYSIS LEGACY (analyze_motifs)

Runs HOMER and kmer enrichment aspects from clip_analysis_legacy to generate useful motif analyses.

# Requirements:
- HOMER
- EMBOSS
- bedtools
- python=2.7
- matplotlib
- seaborn
- numpy
- pandas
- bedtools=2.26
- pybedtools
- htseq
- bx-python
- scikit-learn
- cairosvg=1.0.22

(or run the create_environment.sh script - don't forget to ```source activate analyze_motifs```)

# Installation:
```python setup.py install```

# Usage:
```
analyze_motifs \
--clusters inputs/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.Fc3Pv3.bed \  # BED file (from input normalization or IDR)
--species hg19_v19 \  # sets prefix for finding region bedfiles (see clip_analysis_legacy/data/regions)
--k 6 \  # k for kmer enrichment score search
--nrand 1 \  # sets the number of random bedfile regions to generate to form the background
--out_file outputs/204_01_RBFOX2.svg \  # output motif image file
--pickle outputs/204_01_RBFOX2.pickle \  # output pickle
--genome_fasta inputs/all.fa \  # fasta file
--homer_out outputs/204_01_RBFOX2_homer_out  # directory by which this script uses to store all outputs (will tar this directory at the end)
```

(or see examples/run_analyze_motifs.204_01.sh script)

#### Regarding output parameters...
In the above example, you do NOT need the 204_01_RBFOX2_homer_out/ to exist, but you DO need to make sure outputs/ exists

# Notes about data/regions:

This package expects the following files inside data/regions:

- ${SPECIES}_CDS.bed
- ${SPECIES}_proxintron500.bed
- ${SPECIES}_distintron500.bed
- ${SPECIES}_five_prime_utrs.bed
- ${SPECIES}_three_prime_utrs.bed

where $(SPECIES) is something like: 'hg19' or 'hg19_v19'. Analyze_motifs needs these
 regions to define where to select the random backgrounds from.

### The following species are supported:
- hg19
- hg19_v19
- ce10
- GRCh38
- GRCh38_v24
- mm9

### Using nonsupported species:

If you have the proper GTF file, you can use the create_region_bedfiles
 script inside the annotator package to generate these files, and put
 them into data/regions if they don't exist there already. You will need to
 reinstall the package after doing so.

https://github.com/byee4/annotator/