
GGSOFT/README.md

GGSOFT (Golden Gate DNA Assembly Size-specified Overhang Finding Tool)

Given a user input of overhang size, determine optimal DNA overhang sequences using
a scoring function to maximize pairwise distances between all pairwise combinations
of overhangs, 4**n, where n is the length of the overhang.

Given an overhang size of 4, there are 4**4 = 256 possible unique overhangs.
Given an overhang size of 5, there are 4**5 = 1024 possible unique overhangs.

A scoring table is precalculated given the overhang length, which is used to 
accelerate the scoring calculation by storing pairwise overhand scores in memory.

May later include functionality for gathering sequence homologs to optimize overhang 
locations.
