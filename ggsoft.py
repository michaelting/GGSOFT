#!/usr/bin/python

#==============================================================================
# ggsoft.py
# Michael Ting
# 28 May 2013
#
# Finds top-scoring overhangs for type IIs restriction enzyme digestion 
#   of a user-specified fragment size
#
# Ideas: 
#   - input an overhang size, which calculates pairwise overhang scores.
#       - more similar overhangs receive a worse score, since we want to keep
#         overhangs as different as possible to optimize DNA assembly
#         conditions (especially one-pot reactions)
#       - 4-base overhang = 256 combinations (4^4), n positions ^ 4 bases
#       - use the overhang size to build a table initially to store pairwise
#         score for hashing
#   - input a minimum/maximum fragment size to specify window locations
#
#==============================================================================

from optparse import OptionParser

from Bio import SeqIO

import itertools

# Generator to extract all sequences from FASTA file
# DOES NOT remove the ">" in the name when processed
def process(infile):
    """
    Extracts sequences from a FASTA file.
    
    Input:
        infile  - A FASTA file with a nucleotide sequence
    Output:
        A generator to yield all nucleotide sequences in the file
        
    Note:
        - For the purposes of this program, we want to restrict the number of
          sequences in the input file to ONE (1) sequence. 
        - In the future, we may want to incorporate homology as information to
          find the best locations across a family of gene sequences. This may
          require parsing a FASTA file with multiple sequences.
    """
    name, seq = None, []
    for line in infile:
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line.strip(), []
        else:
            seq.append(line.strip())
    if name:
        yield (name, ''.join(seq))    


def process2(infile):

    record = SeqIO.read(infile, "fasta")
    
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
        A list of lists of indices corresponding to valid overhang combinations
    
    """

    # associate indices with overhang-sized fragments
    subdict = getsubstrings(seq, OHsize)

    # start from a reference index, go from minsize to maxsize

    for index in range(len(subdict.keys())):
        print "hello"

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
        scoretable  - double dictionary table, which can be accessed
                      using scoretable['firstseq']['secondseq']
    Output:
        The sum total score of all pairwise comparisons between overhangs in
        the input list. The lists will later be sorted by score, from which
        we can access the top x percent, e.g. top 10%    
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

"""
# Finds fragments of specified size
def ggsize(seq, minsize, maxsize):

    # dictionary to hold indices of linkers
    # keys: linker nucleotide sequence
    # values: index of first nucleotide in linker
    linkers = dict()

    i = 0
    step = 0
    while i < len(seq):
        link = seq[i:i+4]
        # if linker not yet listed in dictionary
        if not linkers[link]:
            # add link with index to dictionary
            linkers[link] = i
            i += minsize
        elif step > maxsize:
            # deal with abnormalities

    for i in range(0,len(seq), minsize):
        link = seq[i:i+4]
        # if linker not yet listed in dictionary
        if not linkers[link]:
            # add link with index to dictionary
            linkers[link] = i
"""

# Finds a specific number of fragments ---------------------------------------

""" Major issue: need to find sequences inbetween without going off the edges
    or hitting index 0 mod fragnum
    
    Re-write this method to use a window size
"""
def ggnum(seq, fragnum):
    
    length = len(seq)
    stepsize = length / int(fragnum)
    steps = length % int(fragnum)
    # total number of fragments = # overhangs - 1
    totalfrags = int(fragnum)-1 # uses range from 0,...,fragnum

    linkdict = dict()

    """
    This iterates using the index - instead, should iterate using the number
    of fragments to be produced - 1 (3 fragments means 2 linkers)
    """
    """
    # find links
    i = stepsize
    while i < length:
        start = i-2
        end = i+2
        link = seq[start:end]
        # if linker is not already in dictionary
        if link not in linkdict.keys():
            linkdict[link] = start
            i += stepsize
        else:
            i += 1
    """
            
    i = stepsize
    for num in range(totalfrags):
        start = i-2
        end = i+2
        link = seq[start:end]
        # if linker is not already in dictionary
        if link not in linkdict.keys():
            linkdict[link] = start
            i += stepsize
        else:
            i += 1

    # extract links from dictionary
    fraglist = []
    for link, index in linkdict.items():
        frag = (link, index)
        fraglist.append(frag)

    return fraglist

# executes when program run from command line --------------------------------
def main():

    # read command-line arguments --------------------------------------------
    parser = OptionParser()
    parser.add_option("-i", "--in", dest="infile", help="input sequence FASTA file")
    parser.add_option("-o", "--out", dest="outfile", help="name of output FASTA file")
    parser.add_option("-m", "--min", dest="minsize", help="minimum fragment size")
    parser.add_option("-n", "--max", dest="maxsize", help="maximum fragment size")
    parser.add_option("-c", "--count", dest="fragcount", help="number of fragments to produce")
    parser.add_option("-s", "--size", dest="OHsize", help="overhang size in bp")

    """also need an option for enzyme type to overhang bp size"""

    (options, args) = parser.parse_args()

    infile = options.infile
    outfile = options.outfile
    minsize = options.minsize
    maxsize = options.maxsize
    fragcount = options.fragcount
    OHsize = options.OHsize

    template = open(infile)
    newfile = open(outfile, 'w')

    # extract the sequence from the FASTA file -------------------------------
    seqlist = []
    for name, seq in process(template):
        print "name, seq: " + name + ", " + seq
        item = (name, seq)
        #print "item" + item
        seqlist.append(item)

    print seqlist
    
    template.close()

    # Ensure FASTA file has only 1 sequence ----------------------------------
    if len(seqlist) > 1:
        raise IOError("Too many sequences! Only one sequence per file allowed!")
    elif len(seqlist) < 1:
        raise IOError("Not enough sequences! Only one sequence per file allowed!")

    # build the scoring table using the given fragment size
    score_table = buildtable(OHsize)

    # Find the fragments of specified size in the sequence -------------------
    info = seqlist[0][0]
    seq = seqlist[0][1]
    print "info: " + info
    print "seq: " + seq
    #fraglist = goldengate(seq, minsize, maxsize)
    fraglist = ggnum(seq, fragcount)
    
    # Print the fragment results ---------------------------------------------
    print info + "\n"
    for frag in fraglist:
        print str(frag) + "\n"

    # need to write to new file!
    #for frag in fraglist:

if __name__ == "__main__":
    main()
