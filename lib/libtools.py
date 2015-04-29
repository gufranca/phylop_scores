#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Gustavo Starvaggi Franca
Program name: libtools.py
Date: 2014-07-25
Last date modified: 2014-08-01
License: GPL

1. What it does:
   Core functions to be used for phylop simulations or whatever data in the
   same format of phylop scores in bed files. The functions get data from input
   features and allowed/not_allowed regions from bed files to create random or
   flanking regions. Then calculate scores.

2. Input:
   None

3. Output:
   None

4. Usage:
   import libtools

"""


import random
import glob
import os
import subprocess
from numpy import mean
from pybedtools import BedTool



def read_features(features_bed):
    """
    Read features in bed file an create a BedTool object.

    Arg1: features_bed -> A bed file of features.
    Returns -> A BedTool object of the bed file.
    """
    
    features = BedTool(features_bed)
    return features


def get_regions(regions_bed):
    """
    Read a bed file containing allowed regions to generate random intervals.
    It is usually a bed file of genes, transcripts and intron positions.
    
    Arg1: allowed_regions_bed -> A BedTool object of the allowed regions file.
    Returns -> A dictionary of the file containing all regions
    (introns, exons, etc.) by gene and transcript.

    """

    dict_regions = {}
    for region in regions_bed:
        if "_" in region.name:
            name_reg = region.name.split("_")
            gene, transcript = name_reg[0], name_reg[1]        
            if (gene, transcript, region.chrom) in dict_regions:
                dict_regions[(gene,
                              transcript,
                              region.chrom)
                            ].append((int(region.start), int(region.end)))
            else:
                dict_regions[(gene,
                              transcript,
                              region.chrom)
                            ] = [(int(region.start), int(region.end))]
        else:
            dict_regions[region.name] = "None"
    return dict_regions


def get_bed_files(dir_name):
    """
    Get a list of phylop bed files from the directory. Files must follow UCSC
    phylop format name. Ex: chrY.phyloP46way.wigFix.bed.
    
    Arg1: dir_name -> The name of directory to search for bed files.
    Returns -> A dictionary associating chromosome names and file names.

    """

    if not os.path.exists(dir_name):
        raise IOError("Could not find '%s' directory." % dir_name)

    bed_files = glob.glob(dir_name + "/*.bed")
    filenames = [ f.split("/")[-1] for f in bed_files ]
    # dict associating chr name with filename. Ex: dict{chrY:chrY.xx.xx.bed,..}
    dict_files = { f.split(".")[0]: dir_name + "/" + f for f in filenames }
    return dict_files


def check_overlap(feature_string, query_string):
    """
    Check overlap between two bed strings.

    Arg1: feature_string -> target string.
    Arg2: query_string -> query string.
    Returns -> True (if has overlap), False (not overlap).

    """

    feat_bed = BedTool(feature_string, from_string=True)
    query_bed = BedTool(query_string, from_string=True)
    feat_query = feat_bed.intersect(query_bed)
    return bool(feat_query)


def run_bedextract(bed_region, bed_file):
    """
    Using 'bedextract' from BEDOPS package to extract a given region from a
    big query BED file. The query BED file usually is the huge file we want to
    search for the input bed_region.
    
    Arg1: bed_region -> bed region in string format. Ex: "chrX\tstart\tend"
    Arg2: bed_file -> A SORTED bed file containing query regions. 
    Returns -> A BedTool object of all regions in query BED overlapping with
    bed_region.

    """
    
    if bed_region == "NA\tNA\tNA":
        query_regions = ''
    else:
        # echo bed_region
        p1 = subprocess.Popen(['echo', '-e', bed_region], stdout=subprocess.PIPE)
        # pipe echo to bedextract
        p2 = subprocess.Popen(['bedextract', bed_file, '-'], stdin=p1.stdout,
                              stdout=subprocess.PIPE).stdout.read()
        p1.stdout.close()
        query_regions = BedTool(p2, from_string=True)
    return query_regions
    

def calculate_mean_score(query_regions):
    """
    Calculates the mean phylop (or other) scores of a range of features.
    It suposedly must have the phylop (or other) score field.
    
    Arg1: query_regions -> It is usually the returned value of run_bedextract(),
    which is a BedTool object of a range of features. This argument must 
    be a BedTool object. Ex: "chrY    10526   10527   id-26   -1.025000"
    Returns -> A float value of the mean scores of the region.

    """

    # in case of not finding the interval in query bed file.
    if query_regions == '':
        mean_score = "NA"
    else:
        scores = []
        for feature in query_regions:
            scores.append(float(feature.score))
        mean_score =  mean(scores)
    return mean_score


