#! /usr/bin/python

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
"""
def test_parsefasta():
    pass
"""

# checks that scoring table is built correctly
def test_scoretable1():
    """
    test_scoretable1    
    
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
    
    TV = 4
    TS = 1
    ID = 0      
    
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

def test_scoretable2():
    """
    test_scoretable2
    
    input size 2 overhang:
    |NN
     NN|
    should produce 16*16 = size 256 table
        AA  AT  AG  AC  TA  TT  TG  TC  GA  GT  GG  GC  CA  CT  CG  CC
    AA  0   4   1   4   4   8   5   8   1   5   2   5   4   8   5   8 
    AT      0
    AG          0
    AC              0
    TA                  0
    TT                      0
    TG                          0
    TC                              0
    GA                                  0
    GT                                      0
    GG                                          0
    GC                                              0
    CA                                                  0
    CT                                                      0
    CG                                                          0
    CC                                                              0
    """
    
    TV = 4
    TS = 1
    ID = 0
    
    size2table = ggsoft.buildtable(2)
    
    eq_(size2table['AA']['AA'], ID+ID, msg="AA-AA table 2 test failed!")
    eq_(size2table['GT']['GT'], ID+ID, msg="GT-GT table 2 test failed!")
    eq_(size2table['AG']['GA'], TS+TS, msg="AG-GA table 2 test failed!")
    eq_(size2table['AT']['CC'], TV+TS, msg="AT-CC table 2 test failed!")
    eq_(size2table['CC']['TG'], TS+TV, msg="CC-TG table 2 test failed!")
    eq_(size2table['GC']['TA'], TV+TV, msg="GC-TA table 2 test failed!")
    eq_(size2table['AA']['TA'], TV+ID, msg="AA-TA table 2 test failed!")
    eq_(size2table['AA']['GA'], TS+ID, msg="AA-GA table 2 test failed!")

def test_scoretable4():
    """
    test_scoretable4
    
    input size 4 overhang:
    |NNNN
     NNNN|
    should produce 256*256 = size 65536 table
        AAAA
    AAAA
    """
    
    TV = 4
    TS = 1
    ID = 0
    
    size4table = ggsoft.buildtable(4)
    
    eq_(size4table['AAAA']['AAAA'], ID+ID+ID+ID, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['CTGA']['CTGA'], ID+ID+ID+ID, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['AATT']['TTAA'], TV+TV+TV+TV, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['GCGC']['ACAC'], TS+ID+TS+ID, msg="AAAA-AAAA table 4 test failed!")
    eq_(size4table['TCGA']['TAAT'], ID+TV+TS+TV, msg="AAAA-AAAA table 4 test failed!")