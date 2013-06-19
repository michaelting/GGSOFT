#!/usr/bin/python

#==============================================================================
# Unit testing file for ggsoft.py using the nose framework
# Copyright 2013 Michael Ting
# https://github.com/michaelting
#
# Tests can be run using "$ nosetests ggsofttest.py"
# 
# Options:
#   -v          : verbose, displays test names
#   -x          : stop the test at the first failure
#   --nocapture : prevent nose from capturing stdout output
#==============================================================================

import ggsoft
from nose.tools import ok_, eq_, raises

def test_parsefasta():
    """
    test_parsefasta
    
    Ensures that the FASTA file is being read correctly using the Biopython
    package SeqIO.
    """
    testfile1 = open("./tests/test1.fasta", "r")
    seq1 = ggsoft.process(testfile1)
    
    eq_(seq1, "AAAATTTT", msg="fasta parser 2 test 1 reading incorrectly!")
    testfile1.close()

@raises(ValueError)        
def test_parsefail1():
    """
    test_parsefail1
    
    Ensures that the ValueError exception is correctly thrown for multiple
    FASTA sequences in the same input file.
    """
    testfile1 = open("./tests/test2.fasta", "r")
    seq1 = ggsoft.process(testfile1)
    testfile1.close()

@raises(ValueError)
def test_parsefail2():
    """
    test_parsefail2
    
    Ensures that the ValueError exception is correctly thrown for an 
    incorrectly formatted input file
    """
    testfile2 = open("./tests/test_parsefail.fasta", "r")
    seq2 = ggsoft.process(testfile2)
    testfile2.close()
    
def test_parsefail3():
    """
    test_parsefail3
    
    Should correctly parse the file
    """
    testfile3 = open("./tests/test1.fasta", "r")
    seq3 = ggsoft.process(testfile3)
    testfile3.close()

def test_scoretable_size1():
    """
    test_scoretable_size1    
    
    Checks that scoring table is built correctly for 1x1 overhangs  
    
    Input size 1 overhang:
    |N
     N|
    should produce 4*4 = size 16 table
        A   T   G   C
    A   0   4   1   4
    T   4   0   4   1
    G   1   4   0   4
    C   4   1   4   0
    """
    
    TV = 8  # transversion
    TS = 1  # transition
    ID = -10  # identical
    
    size1table = ggsoft.buildtable(1)
    
    eq_(size1table['A']['A'], ID, msg="A-A table 1 test failed!")
    eq_(size1table['A']['T'], TV, msg="A-T table 1 test failed!")
    eq_(size1table['A']['G'], TS, msg="A-G table 1 test failed!")
    eq_(size1table['A']['C'], TV, msg="A-C table 1 test failed!")
    eq_(size1table['T']['A'], TV, msg="T-A table 1 test failed!")
    eq_(size1table['T']['T'], ID, msg="T-T table 1 test failed!")
    eq_(size1table['T']['G'], TV, msg="T-G table 1 test failed!")
    eq_(size1table['T']['C'], TS, msg="T-C table 1 test failed!")
    eq_(size1table['G']['A'], TS, msg="G-A table 1 test failed!")
    eq_(size1table['G']['T'], TV, msg="G-T table 1 test failed!")
    eq_(size1table['G']['G'], ID, msg="G-G table 1 test failed!")
    eq_(size1table['G']['C'], TV, msg="G-C table 1 test failed!")
    eq_(size1table['C']['A'], TV, msg="C-A table 1 test failed!")
    eq_(size1table['C']['T'], TS, msg="C-T table 1 test failed!")
    eq_(size1table['C']['G'], TV, msg="C-G table 1 test failed!")
    eq_(size1table['C']['C'], ID, msg="C-C table 1 test failed!")

def test_scoretable_size2():
    """
    test_scoretable_size2
    
    Checks that scoring table is built correctly for 2x2 overhangs      
    
    input size 2 overhang:
    |NN     |NN
     xx|     xx|
    should produce 16*16 = size 256 table
        AA  AT  AG  AC  TA  TT  TG  TC  GA  GT  GG  GC  CA  CT  CG  CC
    AA  0   4   1   4   4   8   5   8   1   5   2   5   4   8   5   8 
    AT  4   0   4   1   8   4   8   5   5   1   5   2   8   4   8   5
    AG  1   4   0   4   5   8   4   8   2   5   1   5   5   8   4   8
    AC  4   1   4   0   8   5   8   4   5   2   5   1   8   5   8   4
    TA  4   8   5   8   0   4   1   4   4   8   5   8   1   5   2   5
    TT  8   4   8   5   4   0   4   1   8   4   8   5   5   1   5   2   
    TG  5   8   4   8   1   4   0   4   5   8   4   8   2   5   1   5
    TC  8   5   8   4   4   1   4   0   8   5   8   4   5   2   5   1
    GA  1   5   2   5   4   8   5   8   0   4   1   4   4   8   5   8
    GT  5   1   5   2   8   4   8   5   4   0   4   1   8   4   8   5
    GG  2   5   1   5   5   8   4   8   1   4   0   4   5   8   4   8
    GC  5   2   5   1   8   5   8   4   4   1   4   0   8   5   8   4
    CA  4   8   5   8   1   5   2   5   4   8   5   8   0   4   1   4
    CT  8   4   8   5   5   1   5   2   8   4   8   5   4   0   4   1
    CG  5   8   4   8   2   5   1   5   5   8   4   8   1   4   0   4
    CC  8   5   8   4   5   2   5   1   8   5   8   4   4   1   4   0
    """
    
    TV = 8      # transversion
    TS = 1      # transition
    ID = -10    # identical
    P2 = -3     # penalty for length 2 identical string
    
    size2table = ggsoft.buildtable(2)
    
    eq_(size2table['AA']['AA'], ID+ID+P2, msg="AA-AA table 2 test failed!")
    eq_(size2table['GT']['GT'], ID+ID+P2, msg="GT-GT table 2 test failed!")
    eq_(size2table['AG']['GA'], TS+TS, msg="AG-GA table 2 test failed!")
    eq_(size2table['AT']['CC'], TV+TS, msg="AT-CC table 2 test failed!")
    eq_(size2table['CC']['TG'], TS+TV, msg="CC-TG table 2 test failed!")
    eq_(size2table['GC']['TA'], TV+TV, msg="GC-TA table 2 test failed!")
    eq_(size2table['AA']['TA'], TV+ID, msg="AA-TA table 2 test failed!")
    eq_(size2table['AA']['GA'], TS+ID, msg="AA-GA table 2 test failed!")

def test_scoretable_size4():
    """
    test_scoretable_size4
    
    Checks that scoring table is built correctly for 4x4 overhangs      
    
    input size 4 overhang:
    |NNNN
     NNNN|
    should produce 256*256 = size 65536 table
            AAAA    AAAT    AAAG    AAAC    ...
    AAAA    0       4       1       4      
    AAAT    4       0       4       1
    AAAG    1       4       0       4
    AAAC    4       1       4       0
    ...
    """
    
    TV = 8      # transversion
    TS = 1      # transition
    ID = -10    # identical
    # ID penalties accumulate as the sum of powers
    P2 = -3         # penalty for length 2 identical string
    P3 = -3-9       # penalty for length 3 identical string
    P4 = -3-9-27    # penalty for length 4 identical string
    
    size4table = ggsoft.buildtable(4)
    
    eq_(size4table['AAAA']['AAAA'], ID+ID+ID+ID+P4, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['CTGA']['CTGA'], ID+ID+ID+ID+P4, msg="CTGA-CTGA table 4 test failed!")
    eq_(size4table['AATT']['TTAA'], TV+TV+TV+TV, msg="AATT-TTAA table 4 test failed!")
    eq_(size4table['GCGC']['ACAC'], TS+ID+TS+ID, msg="GCGC-ACAC table 4 test failed!")
    eq_(size4table['TCGA']['TAAT'], ID+TV+TS+TV, msg="TCGA-TAAT table 4 test failed!")
    eq_(size4table['AGCT']['ACGT'], ID+TV+TV+ID, msg="AGCT-ACGT table 4 test failed!")
    eq_(size4table['TCTT']['TAAT'], ID+TV+TV+ID, msg="TCTT-TAAT table 4 test failed!")
    eq_(size4table['TACT']['AAAG'], TV+ID+TV+TV, msg="TACT-AAAG table 4 test failed!")
    eq_(size4table['GGCA']['TGAC'], TV+ID+TV+TV, msg="GGCA-TGAC table 4 test failed!")
    eq_(size4table['AGCA']['GAAT'], TS+TS+TV+TV, msg="AGCA-GAAT table 4 test failed!")
    eq_(size4table['TTGG']['TTGG'], ID+ID+ID+ID+P4, msg="TTGG-TTGG table 4 test failed!")
    eq_(size4table['GACT']['GACG'], ID+ID+ID+TV+P3, msg="GACT-GACG table 4 test failed!")
    eq_(size4table['TTAT']['GTAC'], TV+ID+ID+TS+P2, msg="TTAT-GTAC table 4 test failed!")
    
def test_seqsubstr1():
    """
    test_seqsubstr1
    
    Checks that the substrings are correctly parsed from an input sequence
    given a specified overhang size.
    
    "ATGCATG" with overhang size 4 should produce a dictionary as follows:
    {} = {  '0':'ATGC',
            '1':'TGCA',
            '2':'GCAT',
            '3':'CATG'}
    """ 
    
    seq = "ATGCATG"
    size = 4
    
    testdict = {0:'ATGC',
                1:'TGCA',
                2:'GCAT',
                3:'CATG'}
    
    newdict = ggsoft.getsubstrings(seq, size)
    
    for i in range(len(seq)-size+1):
        eq_(newdict[i], testdict[i], msg="Overhang equality at index %d failed!" % i)
    
def test_seqsubstr2():
    """
    test_seqsubstr2
    
    Checks that the substrings are correctly parsed from an input sequence
    given a specified overhang size.
    """
    
    seq = "AAAATTTTGGGGCCCC"
    size = 6
    
    testdict = {0:'AAAATT',
                1:'AAATTT',
                2:'AATTTT',
                3:'ATTTTG',
                4:'TTTTGG',
                5:'TTTGGG',
                6:'TTGGGG',
                7:'TGGGGC',
                8:'GGGGCC',
                9:'GGGCCC',
                10:'GGCCCC'}
    
    newdict = ggsoft.getsubstrings(seq,size)
    
    for i in range(len(seq)-size+1):
        eq_(newdict[i], testdict[i], msg="Overhang equality at index %d failed!" % i)

def test_regionfinder():
    """
    test_regionfinder
    
    Checks that regions are properly found in the input sequence for 
    overhang combination and scoring
    """
    
    seq1 = "AAAATTTT"
           #01234567
    OHsize1 = 4
    minsize1 = 2
    maxsize1 = 6
    correct1 = [[2]]
    rlist1 = ggsoft.find_regions(seq1, OHsize1, minsize1, maxsize1)  
    
    eq_(rlist1, correct1, msg="regionfinder test 1 failed!")

    seq2 = "AAAATTTTGGGGCCCCAAAA"
          # 01234567890123456789
          # aaAATTttggGGCcccaAAA current result
          # aaAATTttGGGGccCCAAaa want this
    OHsize2 = 4
    minsize2 = 2
    maxsize2 = 6
    correct2 = [[2,3,4,5],[8,9,10,11],[14]]
    rlist2 = ggsoft.find_regions(seq2, OHsize2, minsize2, maxsize2)

    eq_(rlist2, correct2, msg="regionfinder test 2 failed!")
    
    testfile1 = open("./tests/P42212.fasta", "r")
    seq3 = ggsoft.process(testfile1)
    OHsize3 = 4
    minsize3 = 20
    maxsize3 = 30
    correct3 = [[20, 21, 22, 23, 24, 25, 26, 27, 28, 29], [50, 51, 52, 53, 54, 55, 56, 57, 58, 59], [80, 81, 82, 83, 84, 85, 86, 87, 88, 89], [110, 111, 112, 113, 114, 115, 116, 117, 118, 119], [140, 141, 142, 143, 144, 145, 146, 147, 148, 149], [170, 171, 172, 173, 174, 175, 176, 177, 178, 179], [200, 201, 202, 203, 204, 205, 206, 207, 208, 209], [230, 231, 232]]
    rlist3 = ggsoft.find_regions(seq3, OHsize3, minsize3, maxsize3)

    eq_(rlist3, correct3, msg="regionfinder test 3 failed!")

def test_validdist():
    """
    test_validdist
    
    Checks that boolean function _valid_distance returns correct values for
    given indices.
    """
    
    OHindexL1 = 4
    OHindexU1 = 8
    mindist1 = 3
    maxdist1 = 9
    
    valid1 = ggsoft._valid_distance(OHindexL1, OHindexU1, mindist1, maxdist1)
    ok_(valid1, msg="valid distance test 1 failed!")
    
    OHindexL2 = 4
    OHindexU2 = 8
    mindist2 = 5
    maxdist2 = 7
    
    valid2 = ggsoft._valid_distance(OHindexL2, OHindexU2, mindist2, maxdist2)
    ok_(not valid2, msg="valid distance test 2 failed!")
    
    OHindexL3 = 4
    OHindexU3 = 7
    mindist3 = 4
    maxdist3 = 8
    
    valid3 = ggsoft._valid_distance(OHindexL3, OHindexU3, mindist3, maxdist3)
    ok_(not valid3, msg="valid distance test 3 failed!")
    
    OHindexL4 = 5
    OHindexU4 = 8
    mindist4 = 4
    maxdist4 = 8
    
    valid4 = ggsoft._valid_distance(OHindexL4, OHindexU4, mindist4, maxdist4)
    ok_(not valid4, msg="valid distance test 4 failed!")

def test_combofinder():
    """
    test_combofinder
    
    Checks that all valid combinations are found correctly
    """
    seq = "AAAATTTTGGGGCCCCAAAA"
    OHsize = 4
    minsize = 2
    maxsize = 6
    validcombos = ggsoft.find_combos(seq, OHsize, minsize, maxsize)

    correct = [(3, 8, 14), (3, 9, 14), (4, 8, 14), (4, 9, 14), (5, 8, 14), (5, 9, 14), (5, 11, 14)]

    eq_(validcombos, correct, msg="combofinder test failed!")

def test_scoreall():
    """
    test_scoreall
    
    Checks that combinations are scored correctly
    """
    seq = "AAAATTTTGGGGCCCCAAAA"
    OHsize = 4
    minsize = 2
    maxsize = 6
    validcombos = ggsoft.find_combos(seq, OHsize, minsize, maxsize)
    subdict = ggsoft.getsubstrings(seq, OHsize)
    scored = ggsoft.scoreall(validcombos, OHsize, subdict)    
    correct = [((4, 8, 14), (68, ['TTTT', 'GGGG', 'CCAA'])), ((4, 9, 14), (68, ['TTTT', 'GGGC', 'CCAA'])), ((5, 9, 14), (68, ['TTTG', 'GGGC', 'CCAA'])), ((3, 8, 14), (68, ['ATTT', 'GGGG', 'CCAA'])), ((3, 9, 14), (68, ['ATTT', 'GGGC', 'CCAA'])), ((5, 8, 14), (43, ['TTTG', 'GGGG', 'CCAA'])), ((5, 11, 14), (43, ['TTTG', 'GCCC', 'CCAA']))]    
    
    eq_(correct, scored, msg="scoreall test failed!")

def test_scorecalc_4():
    """
    test_scorecalc_4
    
    Checks that scores for an overhang list are calculated correctly
    
    "agAGATcaGTCTcaTACGaa" with overhang list ['AGAT','GTCT','TACG']
    should calculate as follows:
    AGAT-GTCT: TS+TV+TV+ID = 1+4+4+0 = 9
    AGAT-TACG: TV+TS+TV+TV = 4+1+4+4 = 13
    GTCT-TACG: TV+TV+ID+TV = 4+4+0+4 = 12
    Sum = 9+13+12 = 34
    """
    
    size = 4
    size4table = ggsoft.buildtable(size)    
    
    OHlist1 = ['AGAT','GTCT','TACG']    
    #testscore1 = 34
    testscore1 = 46
    newscore1 = ggsoft.calc_score(OHlist1,size4table)
    
    eq_(newscore1, testscore1, msg="size 4 score calculation 1 failed!")
    
    OHlist2 = ['AAAA','GATC','GGAT']
    #testscore2 = 21
    testscore2 = 7
    newscore2 = ggsoft.calc_score(OHlist2,size4table)
    
    eq_(newscore2, testscore2, msg="size 4 score calculation 2 failed!")

def test_scorecalc_6():
    """
    test_scorecalc_6
    
    Checks that scores for an overhang list are calculated correctly
    """
    
    TV = 8      # transversion
    TS = 1      # transition
    ID = -10    # identical    
    
    size = 6
    size6table = ggsoft.buildtable(size)
    
    # AAAAAA-TTTTTT: TV*6 = 24
    # AAAAAA-GGGGGG: TS*6 = 6
    # TTTTTT-GGGGGG: TV*6 = 24
    # total score = 54
    OHlist1 = ['AAAAAA','TTTTTT','GGGGGG']
    #testscore1 = 54
    testscore1 = TV*6 + TS*6 + TV*6
    newscore1 = ggsoft.calc_score(OHlist1,size6table)
    
    eq_(newscore1, testscore1, msg="size 6 score calculation 1 failed!")
    
    # AAAAAA-TTTTTT: TV*6 = 24
    # AAAAAA-ACACAC: ID*3+TV*3 = 0+12 = 12
    # AAAAAA-TCTAGT: ID*1+TS*1+TV*4 = 0+1+16 = 17
    # TTTTTT-ACACAC: TS*3+TV*3 = 3+12 = 15
    # TTTTTT-TCTAGT: ID*3+TS*1+TV*2 = 0+1+8 = 9
    # ACACAC-TCTAGT: TV+ID+TV+TV+TS+TS = 0+2+12 = 14
    # total score = 91
    OHlist2 = ['AAAAAA','TTTTTT','ACACAC','TCTAGT']
    #testscore2 = 91
    testscore2 = TV*6 + (ID*3+TV*3) + (ID+TS+TV*4) + (TS*3+TV*3) + (ID*3+TS+TV*2) + (ID+TS*2+TV*3)
    newscore2 = ggsoft.calc_score(OHlist2,size6table)
    
    eq_(newscore2, testscore2, msg="size 6 score calculation 2 failed!")