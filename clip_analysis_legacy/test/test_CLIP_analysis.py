#!/usr/env python

import os
import pandas as pd
import pybedtools
import pytest
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from clip_analysis_legacy import CLIP_analysis as ca

def get_bedtool_1():
    """ returns a bedtool containing 1 interval 10nt long. """
    interval = pybedtools.create_interval_from_list(
        ['chr1','10','20','.','0','+']
    )
    return pybedtools.BedTool([interval])

def get_bedtool_2():
    """
    returns a bedtool containing 1 interval 10nt long
    within the CDS of ENSG00000188976.6
    """
    interval = pybedtools.create_interval_from_list(
        ['chr1','881553','881563','.','0','+']
    )
    return pybedtools.BedTool([interval])

def get_bedtool_3():
    """
    returns a bedtool containing 1 interval 2nt long
    within the CDS of ENSG00000188976.6
    """
    interval = pybedtools.create_interval_from_list(
        ['chr1','881553','881555','.','0','+']
    )
    return pybedtools.BedTool([interval])

def get_bedtool_4():
    """
    returns a bedtool containing 2 intervals 2nt and 10nt long
    within the CDS of ENSG00000188976.6
    """
    interval1 = pybedtools.create_interval_from_list(
        ['chr1','881553','881555','.','0','+']
    )
    interval2 = pybedtools.create_interval_from_list(
        ['chr1', '881553', '881565', '.', '0', '+']
    )

    return pybedtools.BedTool([interval1, interval2])

def get_genome_location():
    return "/projects/ps-yeolab4/genomes/hg19/chromosomes/all.fa"

def get_fasta_dir():
    return "/home/bay001/projects/codebase/clip_analysis_legacy/clip_analysis_legacy/data/regions"

def test_assign_to_regions_2():
    tool = get_bedtool_2()
    regions = ['all', 'cds']
    clusters = "test_assign_to_regions_2"
    fasta_dir = "test_assign_to_regions_2"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )

    assert bed_dict['all']['real'][0].start == 881553
    assert bed_dict['all']['real'][0].end == 881563
    assert bed_dict['cds']['real'][0].start == 881553
    assert bed_dict['cds']['real'][0].end == 881563

def test_assign_to_regions_3():
    tool = get_bedtool_3()
    regions = ['all', 'cds']
    clusters = "test_assign_to_regions_3"
    fasta_dir = "test_assign_to_regions_3"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )

def test_make_fasta_files_from_regions():
    tool = get_bedtool_3()
    regions = ['all', 'cds']
    clusters = "test_make_fasta_files_from_regions"
    fasta_dir = "test_make_fasta_files_from_regions"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )
    ca.make_fasta_files_from_regions(
        cluster_regions=bed_dict,
        clusters=clusters,
        fasta_dir=fasta_dir,
        speciesFA=get_genome_location()
    )

def test_calculate_kmer_diff_1():
    tool = get_bedtool_3()
    regions = ['all', 'cds']
    clusters = "test_calculate_kmer_diff_1"
    fasta_dir = "test_calculate_kmer_diff_1"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )
    ca.make_fasta_files_from_regions(
        cluster_regions=bed_dict,
        clusters=clusters,
        fasta_dir=fasta_dir,
        speciesFA=get_genome_location()
    )

    kmer_diff_result = ca.calculate_kmer_diff(
        kmer_list=kmer_list,
        regions=regions,
        clusters=clusters,
        fasta_dir=fasta_dir
    )

def test_calculate_homer_motifs_1():
    """
    Looks at behavior when looking at one peak < 6nt long.
    """
    tool = get_bedtool_3()
    regions = ['all', 'cds']
    clusters = "test_calculate_homer_motifs_1"
    fasta_dir = "test_calculate_homer_motifs_1"
    homerout = "test_calculate_homer_motifs_1"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )
    ca.make_fasta_files_from_regions(
        cluster_regions=bed_dict,
        clusters=clusters,
        fasta_dir=fasta_dir,
        speciesFA=get_genome_location()
    )
    ca.calculate_homer_motifs(
        kmer_list=kmer_list,
        regions=regions,
        clusters=clusters,
        fasta_dir=fasta_dir,
        homerout=homerout
    )

def test_calculate_homer_motifs_2():
    """
    Looks at behavior when looking at one peak > 6nt long.
    """
    tool = get_bedtool_2()
    regions = ['all', 'cds']
    clusters = "test_calculate_homer_motifs_2"
    fasta_dir = "test_calculate_homer_motifs_2"
    homerout = "test_calculate_homer_motifs_2"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )
    ca.make_fasta_files_from_regions(
        cluster_regions=bed_dict,
        clusters=clusters,
        fasta_dir=fasta_dir,
        speciesFA=get_genome_location()
    )
    ca.calculate_homer_motifs(
        kmer_list=kmer_list,
        regions=regions,
        clusters=clusters,
        fasta_dir=fasta_dir,
        homerout=homerout
    )

def test_calculate_homer_motifs_3():
    """
    Looks at behavior when looking at two peaks.
    peak1 = 2nt long
    peak2 = 10nt long
    """
    tool = get_bedtool_4()
    regions = ['all', 'cds']
    clusters = "test_calculate_homer_motifs_3"
    fasta_dir = "test_calculate_homer_motifs_3"
    homerout = "test_calculate_homer_motifs_3"
    species = "hg19"
    kmer_list = [6]

    bed_dict = ca.assign_to_regions(
        tool=tool,
        clusters=clusters,
        assigned_dir=fasta_dir,
        species=species,
        nrand=1
    )
    ca.make_fasta_files_from_regions(
        cluster_regions=bed_dict,
        clusters=clusters,
        fasta_dir=fasta_dir,
        speciesFA=get_genome_location()
    )
    ca.calculate_homer_motifs(
        kmer_list=kmer_list,
        regions=regions,
        clusters=clusters,
        fasta_dir=fasta_dir,
        homerout=homerout
    )
