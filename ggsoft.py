#!/usr/bin/python

#==============================================================================
# GGSOFT: Golden Gate DNA Assembly Size-specified Overhang Finding Tool
# Copyright 2013 Michael Ting
# https://github.com/michaelting
# Created 28 May 2013
# v1.0 updated 5 June 2013
#
# Finds top-scoring overhangs for type IIs restriction enzyme digestion 
#   of a user-specified fragment size. Overhangs are scored by maximal distance
#   between all overhang pairs using a scoring function which currently
#   uses base substitution (transversion, transition, identical) to correspond
#   to distance, with TV:4, TS:1, ID:0. 
#
# Notes:
#   - Current version uses an oversimplified scoring function. Future versions
#     will explore different scoring functions for higher resolution between
#     scored overhang combinations.
#   - For small fragment sizes, computation of overhang window regions may
#     use up memory and cause the program to fail. Further optimization 
#     may utilize heuristics to reduce memory usage.
#
#==============================================================================

import math, itertools
from argparse import ArgumentParser
from Bio import SeqIO

# Uses the Biopython package SeqIO to extract information from a FASTA file
def process(infile):
    """
    Extracts sequences from a FASTA file using the Biopython package
    SeqIO.
    
    Input:
        infile  - A FASTA file with a single nucleotide sequence.
    Output:
        A string corresponding to the nucleotide sequence in the file.
    Note:
        - For the purposes of this program, we want to restrict the number of
          sequences in the input file to ONE (1) sequence. 
        - In the future, we may want to incorporate homology as information to
          find the best locations across a family of gene sequences. This may
          require parsing a FASTA file with multiple sequences.
    """
    try:
        record = SeqIO.read(infile, "fasta")
    except ValueError:
        raise ValueError("Invalid file format! Please use a FASTA file with one sequence.")
    
    return str(record.seq)

# Construct the scoring table ------------------------------------------------
def buildtable(size):
    """
    Constructs the scoring table of size x size using single-base scores in
    score_dict. Utilizes the itertools.product module to calculate all possible
    nucleotide sequences NNN... of length "size".
    
    Input:
        size    - The length of the overhang in bp
    Output:
        size_score_table    - A matrix of pairwise scores for overhangs of
                              length "size".
    Notes:
        This method is very slow: O(4^n). As we increase the overhang length to
        be greater than 4, the computation becomes extremely slow - 49s for
        a size 6 overhang. We may want to look into rewriting the code to
        speed up this computation, since the matrix is a mirror-image along
        the diagonal.

    input size 1 overhang:
    |N
     N|
    should produce 4*4 = size 16 table
        A   T   G   C
    A   0   4   1   4
    T   4   0   4   1
    G   1   4   0   4
    C   4   1   4   0    
    
    """
    
    TV = 4  # transversion
    TS = 1  # transition
    ID = 0  # identical
    
    # pairwise scores for bases in the same position
    # 1st base --> 2nd base; 1st base is original
    score_dict = {'AA':ID,
                  'AT':TV,
                  'AG':TS,
                  'AC':TV,
                  'TA':TV,
                  'TT':ID,
                  'TG':TV,
                  'TC':TS,
                  'GA':TS,
                  'GT':TV,
                  'GG':ID,
                  'GC':TV,
                  'CA':TV,
                  'CT':TS,
                  'CG':TV,
                  'CC':ID}
    
    # build size n-base scoring table ----------------------------------------
    bases = ['A','T','G','C']    
    
    # find all possible sequences of length "size"
    # size 4 corresponds to "NNNN" = "AAAA","AAAT","AAAG",...,"CCCC"
    overhangs = []
    for seqlst in list(itertools.product(bases, repeat=size)):
        final = ''
        # using itertools.product, seqlst looks like 
        # [['A','A','A','A',],['A','A','A','T'],['A','A','A','G'],...] 
        # so we want to convert each sublist into a string like
        # ['AAAA','AAAT','AAAG',...]
        for letter in seqlst:
            final += letter
        overhangs.append(final)
        
    # construct the overhang table--------------------------------------------
    # --> could optimize this to reduce the number of computations by half by
    # copying the table over in a mirror-image        
        
    # first overhang sequence for comparison
    size_score_table = {}
    for first in overhangs:
        # second overhang sequence to be compared with
        subtable = {}
        for second in overhangs:
            # sum the scores from pairwise comparisons of bases
            # use one-base scoring table to calculate sequence scores
            pairscore = 0
            # go base by base down the overhangs, which should be the same size
            for index in range(size):
                pairstring = first[index] + second[index]
                pairscore += score_dict[pairstring]
            subtable[second] = pairscore
        size_score_table[first] = subtable

    return size_score_table

# Compute all overhang-sized substrings of an input sequence -----------------
def getsubstrings(seq, size):
    """
    Computes all overhang-sized substrings of an input sequence.
    
    Input:
        seq     - The input nucleotide sequence for which we want to find
                  overhangs to use for DNA assembly
        size    - The length in bp of the overhang
    Output:
        A dictionary mapping sequence index --> overhang substring
        
    seq = "ATGCGTA" with overhang size 4 would yield:
    dict = {0:'ATGC',
            1:'TGCG',
            2:'GCGT',
            3:'CGTA'}
    """
    
    subdict = {}
    # slide a window of length "size" across the sequence and record the 
    # substring along with its corresponding initial index
    for pos in range(len(seq)):
        
        start = pos
        end = pos+size
        
        substr = seq[start:end]        
        subdict[pos] = substr
        
        # end of substring sliding window has reached end of sequence
        if end == len(seq):
            break

    return subdict

def find_regions(seq, OHsize, minsize, maxsize):
    """
    Loops through the entire sequence to find regions and store them into
    a list of lists of valid indices for overhang positions.
    
    Input:
        seq     - DNA sequence of interest
        OHsize  - Length of the overhang in bp
        minsize - Minimum fragment size in bp
        maxsize - Maximum fragment size in bp
    Output:
        regionlist - A list of lists of valid indices for overhang positions.    
    """
    
    regionlist = []
    # start from minimum fragment size
    for i in range(minsize,len(seq),maxsize):
        
        region = []
        
        for j in range(maxsize-minsize):
            
            if i >= len(seq)-OHsize-1:
                break
            else:            
                region.append(i)
            i += 1

        # the list "region" is not empty
        if region:
            regionlist.append(region)
        
    return regionlist

# find all valid overhang combinations ---------------------------------------
def find_combos(seq, OHsize, minsize, maxsize):
    """
    Use sequence indices to pull overhang substrings from substrs and produce
    a list of lists of valid overhang combinations

    Input:
        seq     - nucleotide sequence of the input
        OHsize  - length of the overhang in bp
        minsize - minimum size of a DNA fragment
        maxsize - maximum size of a DNA fragment
    Output:
        checked - A list of lists of indices corresponding to valid overhang
                  combinations with their score
    """

    # divide into regions
    regionlist = find_regions(seq, OHsize, minsize, maxsize)  
    
    # find all combinations
    combolist = list(itertools.product(*regionlist))

    checked = []    # holds valid overhang combinations
    # [[0,1],[4,5],[8,9]] becomes
    # [[0,4,8],[0,4,9],[0,5,8],[0,5,9],[1,4,8],[1,4,9],[1,5,8],[1,5,9]]
    for combo in combolist:
        # [0,4,8]
        # checks that 0-->4 and 4-->8 are valid distances; if so, add the 
        # combo [0,4,8] to the list "checked"
        keep = True
        for index in range(len(combo)-1):
            valid = _valid_distance(combo[index], combo[index+1], minsize, maxsize)
            if not valid:
                keep = False
        if keep:
            checked.append(combo)
        
    # return groups of valid indices for overhangs
    return checked

# Checks that regions are within a valid distance from each other ------------
# Alternatively, checks that fragment sizes are within user specifications
def _valid_distance(OHindex1, OHindex2, mindist, maxdist):
    """
    Helper function to determine whether the distance between two overhang 
    window regions satisfies the constraints set by fragment sizes
    specified by the user.
    
    Input:
        OHindex1    - Index of first base in first overhang to compare
        OHindex2    - Index of first base in second overhang to compare
        mindist     - Minimum distance between regions, equivalent to
                      the minimum fragment size
        maxdist     - Maximum distance between regions, equivalent to
                      the maximum fragment size
    Output:
        isValid     - A boolean indicating whether the indices of two overhangs
                      are valid in the sense that they create fragments within
                      the range of minimum to maximum fragment size.
    """
    dist = math.fabs(OHindex2-OHindex1)
    
    isValid = False    
    
    # inclusive, meaning [mindist, maxdist]
    if dist in range(mindist, maxdist+1):
        isValid = True
        
    return isValid

def scoreall(checked, OHsize, subdict):
    """
    Returns a list of valid indices with their corresponding overhangs and
    score, sorted by highest score to lowest score
    
    Input:
        checked - A checked/purged list of overhang index combinations that
                  has been verified to be within a valid distance
        OHsize  - Length of the overhangs in bp
        subdict - A dictionary of (k,v) pairs with keys corresponding to 
                  indices of bases in the sequence and values corresponding to
                  overhang substrings of length OHsize
    Output:
        sorteddict - A dictionary of (k,v) pairs with keys corresponding to
                     index sets and values corresponding to tuples of the
                     overhang set score and the bases in the overhangs 
                     established by the indices.
                     {(206, 413, 620): (17, ['GTTA', 'GTAA', 'CTCA'])}
    """
       
    # create scoring table for overhangs
    scoretable = buildtable(OHsize)
    
    OHcombodict = {}
    # combo looks like [0,5,9]
    for combo in checked:
        overhangs = []
        for seqindex in combo:
            OHseq = subdict[seqindex]
            overhangs.append(OHseq)
        # 'overhangs' has sequences corresponding to indices in 'combo', so
        # [0,5,9] translates to ['AAAA','TTTT','GGGG'], where starting at
        # index 0, a sequence of four (4) bases long is 'AAAA'.
        score = calc_score(overhangs, scoretable)
        value_pair = (score, overhangs)     
        # Key:      (206,413,620)
        # Value:    (17, ['GTTA', 'GTAA', 'CTCA'])
        OHcombodict[combo] = value_pair

    # sort the k,v pairs by score from highest to lowest
    sorteddict = sorted(OHcombodict.items(), key=lambda (k,v): v, reverse=True)
    
    return sorteddict
    
# Find the total score of one overhang combination ---------------------------
def calc_score(OHlist, scoretable):
    """
    Given an overhang list, ['AAAA','GATC','GGAT'], score the list by summing
    all pairwise scores, so:
    
    SUM:
        score('AAAA','GATC') = TS+ID+TV+TV = 1+0+4+4 = 9
        score('AAAA','GGAT') = TS+TS+ID+TV = 1+1+0+4 = 6
        score('GATC','GGAT') = ID+TS+TV+TS = 0+1+4+1 = 6
    = 9+6+9 = 21
    
    For 3 overhangs, we need 2+1 = 3 scores to sum
    For 4 overhangs, we need 3+2+1 = 6 scores to sum
    For 5 overhangs, we need 4+3+2+1 = 10 scores to sum
    For n overhangs, we need 1+2+3+...+n-1 = n*(n-1)/2 = 1/2*(n^2-n) scores

    Input:
        OHlist      - list of overhangs, ['AAAA','GATC','GGAT']
        scoretable  - double dictionary table that contains the pairwise
                      scores for all possible overhang pairs, which can be
                      accessed using scoretable['firstseq']['secondseq']
    Output:
        totalscore  - The sum total score of all pairwise comparisons between 
                      overhangs in the input list.
    """
    
    # Only need to do the upper triangle half of the matrix, otherwise we will
    # be overcalculating the score
    totalscore = 0
    
    for first in range(len(OHlist)):
        
        startsecond = first+1
        seq1 = OHlist[first]
        
        # upper half of matrix only
        for second in range(startsecond,len(OHlist)):
            
            seq2 = OHlist[second]
            score = scoretable[seq1][seq2]
            
            totalscore += score
    
    return totalscore

# executes when program run from command line --------------------------------
def main():

    # read command-line arguments
    parser = ArgumentParser(description="Set parameters for overhang calculation")
    
    parser.add_argument("infile", help="input sequence FASTA file")
    parser.add_argument("outfile", help="name of output file for scored overhangs")
    parser.add_argument("minsize", help="minimum fragment size in bp")
    parser.add_argument("maxsize", help="maximum fragment size in bp")
    parser.add_argument("OHsize", help="overhang size in bp")

    args = parser.parse_args()
    
    infile = args.infile
    outfile = args.outfile
    minsize = int(args.minsize)
    maxsize = int(args.maxsize)
    OHsize = int(args.OHsize)

    template = open(infile)
    newfile = open(outfile, 'w')
    
    # parse the sequence from a FASTA file
    seq = process(template)    
    template.close()
    
    # Find all valid overhang combinations
    combos = find_combos(seq, OHsize, minsize, maxsize)
    
    # Compute index-->overhang dictionary for fast access
    subdict = getsubstrings(seq, OHsize)
    
    # Score and sort all overhang combinations
    scored_list = scoreall(combos, OHsize, subdict)

    # Write score-sorted overhang combinations to a new file
    for scored_combo in scored_list:
        string_combo = str(scored_combo)
        string_combo = string_combo.replace("'","")
        # write to file
        newfile.write("%s\n" % string_combo)

    newfile.close()

if __name__ == "__main__":
    main()