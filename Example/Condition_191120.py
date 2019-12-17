from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 191120 예시 '''
nonTemplateStrand = DNA('TTAGTTACGAGTGGTGGCTGCCCATTGTA').complementReverse()


''' mutation 1 '''
mutationFrom1 = 'w'
mutationTo1 = 'x'
mutation1 = [
    (insert, DNA('GG'))
]
def condition1(dna, poly):
    return poly.length() == 8 and poly.peptideKind() == 8

''' mutation 2 '''
mutationFrom2 = 'x'
mutationTo2 = 'y'
mutation2 = [
    (replace, DNA('yy', allSame), DNA('GG'))
]
def condition2(dna, poly):
    return poly.peptideKind() == 7

''' mutation 3 '''
mutationFrom3 = 'y'
mutationTo3 = 'z'
mutation3 = [
    (replace, DNA(2), DNA(2, allSame))
]
def condition3(dna, poly):
    return poly.haveSequence(dna.getBeforeMutated().getPoly())