#!/usr/bin/python

#==============================================================================
#
# ggsofttest.py
#
# Unit testing file for ggsoft.py using the nose framework
#
#==============================================================================

import ggsoft
from nose.tools import ok_, eq_

# checks that sequence parser functions correctly

def test_parsefasta():
    """
    test_parsefasta
    
    Ensures that the FASTA file is being read correctly
    """   
    eq_(True, False, msg="FASTA file parsed incorrectly!")


# checks that scoring table is built correctly
def test_scoretable_size1():
    """
    test_scoretable_size1    
    
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
    
    TV = 4  # transversion
    TS = 1  # transition
    ID = 0  # identical
    
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
    
    input size 2 overhang:
    |NN
     NN|
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
    
    TV = 4  # transversion
    TS = 1  # transition
    ID = 0  # identical
    
    size2table = ggsoft.buildtable(2)
    
    eq_(size2table['AA']['AA'], ID+ID, msg="AA-AA table 2 test failed!")
    eq_(size2table['GT']['GT'], ID+ID, msg="GT-GT table 2 test failed!")
    eq_(size2table['AG']['GA'], TS+TS, msg="AG-GA table 2 test failed!")
    eq_(size2table['AT']['CC'], TV+TS, msg="AT-CC table 2 test failed!")
    eq_(size2table['CC']['TG'], TS+TV, msg="CC-TG table 2 test failed!")
    eq_(size2table['GC']['TA'], TV+TV, msg="GC-TA table 2 test failed!")
    eq_(size2table['AA']['TA'], TV+ID, msg="AA-TA table 2 test failed!")
    eq_(size2table['AA']['GA'], TS+ID, msg="AA-GA table 2 test failed!")

def test_scoretable_size4():
    """
    test_scoretable_size4
    
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
    
    TV = 4  # transversion
    TS = 1  # transition
    ID = 0  # identical
    
    size4table = ggsoft.buildtable(4)
    
    eq_(size4table['AAAA']['AAAA'], ID+ID+ID+ID, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['CTGA']['CTGA'], ID+ID+ID+ID, msg="CTGA-CTGA table 4 test failed!")
    eq_(size4table['AATT']['TTAA'], TV+TV+TV+TV, msg="AATT-TTAA table 4 test failed!")
    eq_(size4table['GCGC']['ACAC'], TS+ID+TS+ID, msg="GCGC-ACAC table 4 test failed!")
    eq_(size4table['TCGA']['TAAT'], ID+TV+TS+TV, msg="TCGA-TAAT table 4 test failed!")
    eq_(size4table['AGCT']['ACGT'], ID+TV+TV+ID, msg="AGCT-ACGT table 4 test failed!")
    eq_(size4table['TCTT']['TAAT'], ID+TV+TV+ID, msg="TCTT-TAAT table 4 test failed!")
    eq_(size4table['TACT']['AAAG'], TV+ID+TV+TV, msg="TACT-AAAG table 4 test failed!")
    eq_(size4table['GGCA']['TGAC'], TV+ID+TV+TV, msg="GGCA-TGAC table 4 test failed!")
    eq_(size4table['AGCA']['GAAT'], TS+TS+TV+TV, msg="AGCA-GAAT table 4 test failed!")
    eq_(size4table['TTGG']['TTGG'], ID+ID+ID+ID, msg="TTGG-TTGG table 4 test failed!")
    
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
        
def test_scorecalc():
    """
    test_scorecalc
    
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
    testscore1 = 34
    newscore1 = ggsoft.calc_score(OHlist1,size4table)
    
    eq_(newscore1, testscore1, msg="size 4 score calculation 1 failed!")
    
    OHlist2 = ['AAAA','GATC','GGAT']
    testscore2 = 21
    newscore2 = ggsoft.calc_score(OHlist2,size4table)
    
    eq_(newscore2, testscore2, msg="size 4 score calculation 2 failed!")
    
def test_scorecalc2():
    """
    test_scorecalc2
    
    Checks that scores for an overhang list are calculated correctly
    """
    
    size = 6
    size6table = ggsoft.buildtable(size)
    
    # AAAAAA-TTTTTT: TV*6 = 24
    # AAAAAA-GGGGGG: TS*6 = 6
    # TTTTTT-GGGGGG: TV*6 = 24
    # total score = 54
    OHlist1 = ['AAAAAA','TTTTTT','GGGGGG']
    testscore1 = 54
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
    testscore2 = 91
    newscore2 = ggsoft.calc_score(OHlist2,size6table)
    
    eq_(newscore2, testscore2, msg="size 6 score calculation 2 failed!")
    
def test_combofinder1():
    """
    test_combofinder1
    
    Checks that all valid combinations are found correctly
    """
    
    # fragment size between 4-8bp    
    seq = "aaaaTTTTaaaa"
    """    
    TTTT
    """
    
    