#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Gustavo Starvaggi Franca
Program name: features.py
Date: 2014-07-22
Last date modified: 2014-08-04
License: GPL

1. What it does:
   It defines a class for BED features, which will be tested for random or
   flanking intervals, getting their scores from query BED files.

2. Input:
   None

3. Output:
   None

4. Usage:
   import features.py

"""


import random
from numpy import mean
from pybedtools import BedTool
from libtools import(read_features,
                     get_regions,
                     calculate_mean_score,
                     run_bedextract,
                     check_overlap,
                     get_bed_files)


class Feature(object):
    """
    Creates an object for a feature entry, which is usually a single line of
    a BedTool object, where we want to do simulate phyloP or whatever scores
    extracted from a big Bed file containing such scores.

    """

    def __init__(self, feature):
        """
        Initialize the attributes from a BedTool entry.

        """

        self.chrom = feature.chrom
        self.start = int(feature.start)
        self.end = int(feature.end)
        self.name = feature.name
        self.strand = feature.strand
        self.fields = feature.fields
        self.size = int(feature.end) - int(feature.start)

    def flanking_regions(self, window_r=1, window_l=1):
        """
        Define the flanking regions around feature start and end.
        
        Returns -> Two strings, containing right and left intervals. If they are
        valid intervals, will be checked in other by flanking_simulation().

        """

        right_interval = (self.end + window_r, self.end + self.size + window_r)
        left_interval = ((self.start - window_l) - self.size, self.start - window_l)
        right = "%s\t%d\t%d" % (self.chrom, right_interval[0],right_interval[1])
        left = "%s\t%d\t%d" % (self.chrom, left_interval[0], left_interval[1])
        return right, left
        
    def flanking_simulation(self, not_allowed_regions_bed, query_bed):
        """
        Performs a simulation on features, selecting flanking intervals based on
        input feature, extracts the region in query_bed containing scores and
        calculates the mean score over the flanking region. It tries to get 
        valid flanking regions (i.e, those that do not overlaps with regions in 
        not_allowed_regions_dict and have scores annotated in query BEDs)
        If the flanking region overlaps to not allowed regions or do not have
        scores, the score will be 'NA'.

        Arg1: not_allowed_regions_bed -> Must be a sorted bed file containing 
        regions to filter out.
        Arg2: query_bed -> A filename of a BED file, containing all regions and
        scores.

        Returns -> A float score, which is the mean score of the flanking region.
        If only right or left borders have scores, the score will be the one of
        them. If both flanking regions have scores, it will be the mean of both,
        and if no flanking regions have scores, it will be NA.

        """

        # get right and left flanking regions
        right_flank, left_flank = self.flanking_regions()
        # check if the flanking region intersects with not allowed regions
        intersect_r = run_bedextract(right_flank, not_allowed_regions_bed)
        intersect_l = run_bedextract(left_flank, not_allowed_regions_bed)
        
        # testing right and left flanking regions and calculate scores
        # intersected regions not empty, means that overlap with not allowed 
        # regions so we do not want a score for that.
        if (intersect_r != '') and (intersect_l != ''):
            score = "NA"
        elif (intersect_r != '') and (intersect_l == ''):
            left_feature = run_bedextract(left_flank, query_bed)
            score = calculate_mean_score(left_feature)
        elif (intersect_r == '') and (intersect_l != ''):
            right_feature = run_bedextract(right_flank, query_bed)
            score = calculate_mean_score(right_feature)
        elif (intersect_r == '') and (intersect_l == ''):
            left_feature = run_bedextract(left_flank, query_bed)
            right_feature = run_bedextract(right_flank, query_bed)
            score_l = calculate_mean_score(left_feature)
            score_r = calculate_mean_score(right_feature)
            if score_l != "NA" and score_r != "NA":
                score = mean([score_l, score_r])
            elif score_l != "NA" and score_r == "NA":
                score = score_l
            elif score_l == "NA" and score_r != "NA":
                score = score_r
            elif score_l == "NA" and score_r == "NA":
                score = "NA"
        return score
        
    def random_regions(self, allowed_regions):
        """
        By using feature start and end, select random intervals of the same
        length from 'allowed_regions'. These regions usually are intronic 
        regions. Before using this method, you have to check if the feature name
        corresponds to the name of allowed regions. Ex. For an intragenic miRNA,
        get Host Gene and Transcript names and get the corresponding allowed
        regions for the same Gene and Transcript. This information is got by
        get_allowed_regions() function, which collects these info in a dict.

        It tries to get valid random regions (i.e, regions that do not overlap
        with the feature itself and are within the real interval. If this fails,
        random region will be "NA NA NA".

        Arg1: allowed_regions -> The value of the dictionary generated by
        allowed_regions(). Ex: [(10, 200), (300, 500)...]
        Returns -> A string of a random region, chosen from these allowed
        intervals. The random region do not surpasses the end of an allowed
        interval and do not intersects with the feature itself.

        """
        
        feature_string = "%s\t%s\t%s" % (self.chrom, self.start, self.end)
        attempts = 0
        MAX_ATTEMPTS = 100
        while True:
            attempts += 1
            regions_copy = allowed_regions[:]
            random.shuffle(regions_copy)
            # just pick the first one
            region_start = regions_copy[0][0]
            region_end = regions_copy[0][1]
            random_start = random.randrange(region_start, region_end)
            random_end = random_start + self.size
            random_string = "%s\t%s\t%s" % (self.chrom, random_start, random_end)
            # check if random interval overlaps to feature itself.
            intersect = check_overlap(feature_string, random_string)
            # check if the random end do not surpasses the allowed end and
            # do not intersect with the feature itself.
            if (random_end < region_end) and (intersect == False):
                random_region = random_string
                break
            elif attempts == MAX_ATTEMPTS:
                random_region = "NA\tNA\tNA"
                break    
        return random_region

    def random_intragenic_simulation(self, allowed_regions_dict, query_bed):
        """
        Performs a simulation on features, selecting random intervals based on
        input feature, extracts the region in query_bed containing scores and
        calculates the mean score over the selected region. It tries to get 
        valid random regions (i.e, those annotated and having scores in query
        BEDs) for 100 times. If the region is empty, the score will be 'NA'.

        Arg1: allowed_regions_dict -> A dictionary containing allowed regions to
        generate random intervals.
        Arg2: query_bed -> A filename of a BED file, containing all regions and
        scores.

        Returns -> A float score, which is the mean score of the random region.

        """
        
        # check intragenic name
        if "_" in self.name:
            names = self.name.split("_")
        else:
            score = "NameError1"
            return score
        # get the combination of Gene, Transcript, Chromosome and Strand in
        # bed regions file dictionary (-b option)
        key = (names[1], names[2], self.chrom)
        if key in allowed_regions_dict:
            allowed_regions = allowed_regions_dict[key]
        else:
            score = "NameError2"
            return score

        attempts = 0
        MAX_ATTEMPTS = 100
        while True:
            attempts += 1
            random_region = self.random_regions(allowed_regions)
            random_features = run_bedextract(random_region, query_bed)
            if random_features != '':
                score = calculate_mean_score(random_features)
                break
            elif attempts == MAX_ATTEMPTS:
                score = calculate_mean_score(random_features)
                break
        return score


    def random_flanking_regions(self, not_allowed_regions_bed,
                                window_r=10000, window_l=10000):
            """
            By using feature start and end, select random intervals of the same
            length from 'allowed_regions'. These regions usually are intronic 
            regions. It tries to get valid random regions (i.e, regions that do 
            not overlap with the feature itself and are within the real 
            interval. The interval to generate random regions is the distance
            between the start/end and the distance of the first overlapping
            feature. If there are no overlapping features, the interval will be
            window_r and window_l for down and upstream flanking regions.
            If this fails, random region will be "NA NA NA".

            Arg1: not_allowed_regions -> A sorted BED file containing regions
            that cannot overlap with the generated random intervals.
            Arg2/3: window_r/l -> Downstream and upstream windows in respect to
            the feature coordinates.
            Returns -> A string of a random region, not overlapping to not
            allowed regions, do not overlapping with the feature itself and do
            not surpasses upstream limits.

            """
            
            feature_string = "%s\t%s\t%s" % (self.chrom, self.start, self.end)
            right = sorted([self.end + 1, self.end + window_r])
            # if upstream window is larger than start, start of left will be 0
            if window_l >= self.start:
                left = sorted([self.start -1, 0])
            else:
                left = sorted([self.start - 1, self.start - window_l])
            
            attempts = 0
            MAX_ATTEMPTS = 100
            while True:
                attempts += 1
                flanking = random.choice([right, left])
                flank = "%s\t%d\t%d" % (self.chrom, flanking[0], flanking[1])
                gene_overlap = run_bedextract(flank, not_allowed_regions_bed)
                # if flanking range overlaps to gene regions, get a shorter 
                # interval.
                if gene_overlap != '':
                    # if the upstream border was chosen, get upstream 
                    # available region.
                    if flanking[1] < self.start:
                        g_index = len(gene_overlap) - 1
                        distance = self.start - int(gene_overlap[g_index].end)
                        random_start = random.randrange(self.start - distance,
                                                        self.start)
                        random_end = random_start + self.size
                    # if the downstream border was chosen, get downstream 
                    # available region
                    elif flanking[1] > self.end:
                        distance = int(gene_overlap[0].start) - self.end
                        random_start = random.randrange(self.end,
                                                        self.end + distance)
                        random_end = random_start + self.size
                # if do not overlap with gene regions, use the whole interval
                else:
                    random_start = random.randrange(flanking[0], flanking[1])
                    random_end = random_start + self.size
                
                random_string = "%s\t%s\t%s" % (self.chrom, random_start, 
                                                random_end)
                # check if random interval overlaps to feature itself
                intersect = check_overlap(feature_string, random_string)
                # check if the random end do not surpasses the allowed end and
                # do not intersect with the feature itself.
                if (random_end < flanking[1]) and (intersect == False):
                    random_region = random_string
                    break
                elif attempts == MAX_ATTEMPTS:
                    random_region = "NA\tNA\tNA"
                    break
            return random_region

    def random_flanking_simulation(self, not_allowed_regions_bed, query_bed,
                                   window_r=10000, window_l=10000):
        """
        Performs a simulation on features, selecting random intervals based on
        input feature, extracts the region in query_bed containing scores and
        calculates the mean score over the selected region. It tries to get 
        valid random regions (i.e, those annotated and having scores in query
        BEDs) for 100 times. If the region is empty, the score will be 'NA'.

        Arg1: allowed_regions_dict -> A dictionary containing allowed regions to
        generate random intervals.
        Arg2: query_bed -> A filename of a BED file, containing all regions and
        scores.

        Returns -> A float score, which is the mean score of the random region.

        """

        attempts = 0
        MAX_ATTEMPTS = 100
        while True:
            attempts += 1
            random_region = self.random_flanking_regions(
                                                        not_allowed_regions_bed,
                                                        window_r, window_l)

            random_features = run_bedextract(random_region, query_bed)
            # Get score for non empty query features        
            if random_features != '':
                score = calculate_mean_score(random_features)
                break
            elif attempts == MAX_ATTEMPTS:
                score = calculate_mean_score(random_features)
                break
        return score


