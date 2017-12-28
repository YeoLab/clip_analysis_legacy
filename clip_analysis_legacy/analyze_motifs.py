#!/usr/bin/env python

__author__ = ['gpratt','byee']

import cPickle as pickle
import tarfile
from argparse import ArgumentParser

import matplotlib as mpl
from matplotlib import rc

import CLIP_analysis_display
from CLIP_analysis import *

mpl.rcParams['svg.fonttype'] = 'none'
rc('text', usetex=False)
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})


def main(bedtool, species, motifs=[], k=[6], nrand=3,
         genome_location=None,
         motif_location=os.getcwd(),
         out_pickle=None, out_file=None, homer_outdir=None
         ):

    """

    Runs all analysies

    one thing to do is make graphs fail gracefully

    """

    print "starting"
    #gets clusters in a bed tools + names species
    clusters = os.path.basename(bedtool)
    species = species
    short_species = species.split("_")[0]

    out_dict = {}
    #In case names aren't unique make them all unique
    clusters_bed = pybedtools.BedTool(make_unique(pybedtools.BedTool(bedtool))).saveas()

    if len(clusters_bed) <= 1:
        print "No Clusters, killing now to save time"
        return

    #makes output file names
    clusters = str.replace(clusters, ".BED", "")


    # k = [int(x) for x in k]

    #sets up output dirs
    make_dir(homer_outdir)

    assigned_dir = os.path.join(homer_outdir, "assigned")
    misc_dir = os.path.join(homer_outdir, "misc")
    fasta_dir = os.path.join(homer_outdir, "fasta")
    homerout_base = os.path.join(homer_outdir, "homer")
    make_dir(homerout_base)
    homerout = os.path.join(homerout_base, clusters)

    make_dir(assigned_dir)
    make_dir(misc_dir)
    make_dir(fasta_dir)
    make_dir(homerout)

    assigned_regions, regions = regions_generator()
    genic_region_sizes = count_genomic_region_sizes(assigned_regions, species)
    cluster_regions = assign_to_regions(tool=clusters_bed, clusters=clusters,
                                        assigned_dir=assigned_dir,
                                        species=species, nrand=nrand)
    region_sizes = get_sizes(cluster_regions)

    if genome_location is not None:
        make_fasta_files_from_regions(cluster_regions, clusters, fasta_dir, genome_location)
        calculate_homer_motifs(k, regions, clusters, fasta_dir, homerout)
        kmer_results = calculate_kmer_diff(k, regions, clusters, fasta_dir)
    else:
        print "No genome fasta file provide, motif identification will not be performed"

    motif_distances = []
    try:
        if motifs:
            motif_distances = generate_motif_distances(
                cluster_regions, region_sizes, motifs, motif_location,
                short_species
            )
    except:
        pass

    out_dict["kmer_results"] = kmer_results
    out_dict["motifs"] = motifs
    out_dict["motif_distances"] = motif_distances
    out_dict['homerout'] = homerout
    out_dict['regions'] = regions
    out_dict["region_sizes"] = region_sizes
    out_dict['genic_region_sizes'] = genic_region_sizes
    with open(out_pickle, 'w') as p:
        pickle.dump(out_dict, file=p)
    print "file saved"

    visualize(out_pickle, out_file)


    try:
        if motifs:
            motif_fig = CLIP_analysis_display.plot_motifs(motif_distances)
            motif_fig.savefig(out_file + '.motif_distribution.svg')
    except:
        pass

    make_tarfile(
        os.path.join(os.path.abspath(os.path.join(homer_outdir, os.pardir)),
                     '{}.tar.gz'.format(os.path.basename(homer_outdir))),
        homer_outdir
    )

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

def visualize(pickle, out_file):
    with open(pickle) as p:
        clip_viz = CLIP_analysis_display.ClipVisualization(p)
    qc_fig = clip_viz.CLIP_QC_figure()
    qc_fig.savefig(out_file)


def call_main():
    parser = ArgumentParser()

    parser.add_argument(
        "--clusters",
        help="BED file of clusters"
    )
    parser.add_argument(
        "--species",
        help = "genome version",
        required=False,
        default='hg19'
    )
    parser.add_argument(
        "--motifs",
        help="Motifs to use (files of motifs give must"
             " exist in motif_directory directory)",
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--k",
        help="k-mer and homer motif ananlysis",
        nargs='+',
        default=[6],
        required=False
    )
    parser.add_argument(
        "--nrand",
        default=3,
        help="selects number of times to randomly sample genome",
        type=int,
        required=False
    )
    parser.add_argument(
        "--pickle",
        help="destination of the pickle output",
        required=True
    )
    parser.add_argument(
        "--out_file",
        help="destination of the image summary output",
        required=True
    )
    parser.add_argument(
        "--homer_out",
        help="destination of the homer intermediate/output files",
        required=True
    )
    ##Below here are critical files that always need to be referenced
    parser.add_argument(
        "--genome_fasta",
        help="location of all.fa file for genome of interest",
        default=None,
        required=True,
    )
    parser.add_argument(
        "--motif_directory",
        dest="motif_location",
        help="directory of pre-computed motifs for analysis",
        default=os.getcwd()
    )

    args = parser.parse_args()

    #error checking
    if args.clusters is None or args.species is None:
        parser.print_help()
        exit()
    main(
        bedtool=args.clusters, species=args.species,
        motifs=args.motifs, k=args.k, nrand=args.nrand,
        genome_location=args.genome_fasta,
        motif_location=args.motif_location,
        out_pickle=args.pickle, out_file=args.out_file,
        homer_outdir=args.homer_out
    )

if __name__ == "__main__":
    call_main()
