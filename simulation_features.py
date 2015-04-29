#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Gustavo Starvaggi Franca
Program name: simulation_features.py
Date: 2014-07-26
Last date modified: 2014-08-04
License: GPL

1. What it does:
   Performs phyloP (or other scoring system) on random or flanking regions of
   input features. A feature object is created and the simulation methods are
   called, depending of with type of simulation we want to perform.

2. Input:
   A Bed file containing features to be simulated. Sorted Bed files of each
   chromosome, containing scores and A sorted Bed file containing regions to
   look for random regions based on the feature file. This file is use to search
   for regions (if feature is intragenic: use therandom_intragenic_simulation()
   method or to  filter out unwanted regions (if feature is intergenic: use
   flanking_simulation() method).

3. Output:
   Feature name and mean scores for flanking or random simulations.

4. Usage:
   python simulation_features.py --help

"""


import sys
import argparse
import lib
from lib.features import Feature
from lib.libtools import read_features, get_bed_files, get_regions


def call_flanking_simulation(features, bed_files, not_allowed_regions_bed):
    """
    Perform simulations for all features, calling flanking_simulation().
    
    Arg1: features -> BedTool object for all features.
    Arg2: bed_files -> dictionary of chromosome and query BED file names.
    Arg3: not_allowed_regions_dict -> dictionary of not allowed regions, which 
    is the search space to avoid flanking regions.

    Returns -> None. Just prints out the output.

    """
 
    for f in features:
        feature = Feature(f) # create feature object
        # get the correct bed file for specific chromosome
        try:
            query_bed = bed_files[feature.chrom]
        except KeyError:
            print "Could not find a BED file for %s." % (feature.chrom)
            continue
        
        score = feature.flanking_simulation(not_allowed_regions_bed, query_bed)
        out = "%s\t%s" % (feature.name, str(score))
        print out


def call_random_intragenic_simulation(features, bed_files,
                                      allowed_regions_dict, number=1):
    """
    Perform simulations for all features, calling random_simulation_intragenic()
    
    Arg1: features -> BedTool object for all features.
    Arg2: bed_files -> dictionary of chromosome and query BED file names.
    Arg3: allowed_regions_dict -> dictionary of allowed regions, which is the
    search space to generate random intervals and get scores.
    Arg4: number -> Number of simulations to be performed.

    Returns -> None. Just prints out the output.

    """
    
    for f in features:
        feature = Feature(f) # create feature object
        # get the correct bed file for specific chromosome
        try:
            query_bed = bed_files[feature.chrom]
        except KeyError:
            print "Could not find a BED file for %s." % (feature.chrom)
            continue

        scores = []
        for n in range(0, number):
            score = feature.random_intragenic_simulation(allowed_regions_dict,
                                                         query_bed)
            scores.append(score)
        # prepare output
        s = "\t".join([ str(i) for i in scores ])
        out = "%s\t%s" % (feature.name, s)
        print out


def call_random_flanking_simulation(features, bed_files, 
                                    not_allowed_regions_bed, number=1,
                                    window_r=10000, window_l=10000):
    """
    Perform simulations for all features, calling random_flanking_simulation()
    
    Arg1: features -> BedTool object for all features.
    Arg2: bed_files -> dictionary of chromosome and query BED file names.
    Arg3: not_allowed_regions_dict -> dictionary of not allowed regions, which 
    is the search space to avoid flanking regions.
    Arg4: number -> Number of simulations to be performed.
    Arg5/6: window_r/l -> Downstream and upstream windows in respect to the
    feature coordinates.
    Returns -> None. Just prints out the output.

    """

    for f in features:
        feature = Feature(f) # create feature object
        # get the correct bed file for specific chromosome
        try:
            query_bed = bed_files[feature.chrom]
        except KeyError:
            print "Could not find a BED file for %s." % (feature.chrom)
            continue
        
        scores = []
        for n in range(0, number):
            score = feature.random_flanking_simulation(not_allowed_regions_bed,
                                                       query_bed,
                                                       window_r,
                                                       window_l)
            scores.append(score)
        # prepare output
        s = "\t".join([ str(i) for i in scores ])
        out = "%s\t%s" % (feature.name, s)
        print out


def main():
    """
    Get arguments and call functions to perform score simulations on features.

    """

    parser = argparse.ArgumentParser(description="""Performs simulations on 
            phyloP scores (or whatever score) based on  input BED coordinates 
            and querying a BED file containing scores for 
            regions. It is possible to perform two types of simulations.\n

            a) Selecting random intervals based on the input, limiting the 
            search space based on another BED files containing the allowed 
            regions. These random intervals are used to retrieve the
            corresponding Scores in the query BED file. The output is the mean
            score for the selected interval, which can be repeated -n times.\n
            
            b) Selecting the flanking regions based on the input coordinates.
            The mean Scores of the flanking regions are then calculated. It is
            possible to pass regions NOT allowed to overlap with the flanking
            region. In this case, the score cannot be computed for a particular
            feature.
            
            !Warning: Call this program in an external bash loop to avoid
            'Too many files open' error.""") 

    parser.add_argument("-i", dest="features_bed", required=True,
                        help=""""Input file, with features to be simulated.
                        If features are intragenic, fields must be separated
                        by '_'. Ex: feat1_ENSGX_ENSTX... If intragenic names
                        does not follow this pattern scores will be 
                        'NameError1'. This only matters for [-r] option.""")
    parser.add_argument("-d", dest="dirname_bed", required=True,
                        help="""Name of the directory where the BED files per
                        chromosome containing scores are stored. File pattern
                        must be 'chrXX.[...].bed'. *** FILES MUST BE SORTED.""")
    parser.add_argument("-b", dest="regions_bed", required=True,
                        help="""BED file with regions to be considered for
                        searching [-r] or to be filtered out [-f | -rf]. If the 
                        file is gene features (e.g exons or introns), fields 
                        must be separated by '_'. Ex: ENSGX_ENSTX... If names 
                        does not follow this pattern or for some reason, the
                        feature names (gene, transcript, chrom, strand) could
                        be find in this file, scores will be 'NameError2'. 
                        This issue only matters for [-r] option.
                        *** THIS FILE MUST BE SORTED.""")
    
    simulation_group = parser.add_mutually_exclusive_group(required=True)
    simulation_group.add_argument("-f", "--flanking", dest="flanking",
                                  action='store_true', 
                                  help="Flanking simulation.")
    simulation_group.add_argument("-r", "--random", dest="random",
                              action="store_true", help="""Random simulation, 
                              usually  for intragenic regions.""")
    simulation_group.add_argument("-rf", "--random_flank", dest="random_flank",
                              action="store_true", help="""Random flanking 
                              simulation, usually for intergenic regions.""")

    parser.add_argument("-n", "--number", dest="number", type=int, default=1,
                        help="""Number of simulations to be performed. Default
                        = 1. Only accepted with -r or -rf options.""")
    parser.add_argument("-wr", "--down_window", dest="window_down", type=int,
                        default=10000, help="""Downstream window size.
                        Default = 10000. Only accepted with -r or -rf option.""")
    parser.add_argument("-wl", "--up_window", dest="window_up", type=int, 
                        default=10000, help="""Upstream window size. 
                        Default = 10000. Only accepted with -r or -rf option.""")
    
    args = parser.parse_args()

    # checking options
    if args.flanking == True and args.number != 1:
        parser.error("-n is only accepted with -r or -rf options.")
    if (args.flanking == True or args.random == True) and \
       (args.window_down != 10000 or args.window_up != 10000):
        parser.error("-wr or -wl are only accepted with -rf option.")

    #-------------------------------#
    # Call functions and get output #
    #-------------------------------#

    # get features to be tested
    features = read_features(args.features_bed)
    bed_files = get_bed_files(args.dirname_bed)

    # Flanking simulations
    if args.flanking:
        not_allowed_regions_bed = args.regions_bed
        call_flanking_simulation(features, bed_files, not_allowed_regions_bed)
 
    # Random simulations
    elif args.random:
        allowed_regions_bed = read_features(args.regions_bed)
        allowed_regions_dict = get_regions(allowed_regions_bed)
        call_random_intragenic_simulation(features, bed_files,
                                          allowed_regions_dict,
                                          number=args.number)
    # Random flanking simulations
    elif args.random_flank:
        not_allowed_regions_bed = args.regions_bed
        call_random_flanking_simulation(features, bed_files,
                                        not_allowed_regions_bed,
                                        number=args.number,
                                        window_r=args.window_down,
                                        window_l=args.window_up)
        

if __name__ == "__main__":
    main()


#!TODO HOW TO DEAL WITH MULTIPLE CHROMOSOMES IN A SINGLE BED FILE? THIS MUST
# be fixed for small analysis...

