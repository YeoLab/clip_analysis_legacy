#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='analyze_motifs',
    version='0.0.1',
    packages=['clip_analysis_legacy'],
    package_dir={
        'clip_analysis_legacy':'clip_analysis_legacy'
    },
    package_data = {
        'clip_analysis_legacy': [
            'data/regions/*.bed'
        ]
    },
    include_package_data=True,
    url='',
    license='',
    author=['gpratt','byee4'],
    author_email='',
    description='original clip analysis, now just runs homer and kmer enrichment.',
    entry_points = {
        'console_scripts': [
            'analyze_motifs = clip_analysis_legacy.analyze_motifs:call_main',
        ]
    }
)
