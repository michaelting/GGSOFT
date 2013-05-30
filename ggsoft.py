#! /usr/bin/python

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

import sys
import re
from optparse import OptionParser

# Generator to extract all sequences from FASTA file
# DOES NOT remove the ">" in the name when processed
def process(infile):
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

# Finds a specific number of fragments

""" Major issue: need to find sequences inbetween without going off the edges
    or hitting index 0 mod fragnum
"""
def ggnum(seq, fragnum):
    
    length = len(seq)
    stepsize = length / int(fragnum)
    steps = length % int(fragnum)
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

def main():

    parser = OptionParser()
    parser.add_option("-i", "--in", dest="infile", help="input sequence FASTA file")
    parser.add_option("-o", "--out", dest="outfile", help="name of output FASTA file")
    parser.add_option("-m", "--min", dest="minsize", help="minimum fragment size")
    parser.add_option("-n", "--max", dest="maxsize", help="maximum fragment size")
    parser.add_option("-c", "--count", dest="fragcount", help="number of fragments to produce")

    """also need an option for enzyme type to overhang bp size"""

    (options, args) = parser.parse_args()

    infile = options.infile
    outfile = options.outfile
    minsize = options.minsize
    maxsize = options.maxsize
    fragcount = options.fragcount

    template = open(infile)
    newfile = open(outfile, 'w')

    # extract the sequence from the FASTA file
    seqlist = []
    for name, seq in process(template):
        print "name, seq: " + name + ", " + seq
        item = (name, seq)
        #print "item" + item
        seqlist.append(item)

    print seqlist

    # Ensure FASTA file has only 1 sequence
    if len(seqlist) > 1:
        raise IOError("Too many sequences! Only one sequence per file allowed!")
    elif len(seqlist) < 1:
        raise IOError("Not enough sequences! Only one sequence per file allowed!")

    # Find the fragments of specified size in the sequence
    info = seqlist[0][0]
    seq = seqlist[0][1]
    print "info: " + info
    print "seq: " + seq
    #fraglist = goldengate(seq, minsize, maxsize)
    fraglist = ggnum(seq, fragcount)
    
    # Print the fragment results
    print info + "\n"
    for frag in fraglist:
        print str(frag) + "\n"

    # need to write to new file!
    #for frag in fraglist:

if __name__ == "__main__":
    main()
