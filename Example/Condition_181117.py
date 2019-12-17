from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 181117 예시 '''
nonTemplateStrand = DNA('GCATGTTACTCAGCGCTCGCAACTAGCATACATGT'), DNA('GCATGTTACTCAGCGCTCGCAACTAGCATACATGT').complementReverse()


''' mutation 1 '''
mutationFrom1 = 'w'
mutationTo1 = 'x'
mutation1 = [
    (delete, DNA(1), range(8, 11))
]
def condition1(dna, poly):
    return poly.havePeptide(LEU)

''' mutation 2 '''
mutationFrom2 = 'w'
mutationTo2 = 'y'
mutation2 = [
    (replace, DNA(1))
]
def condition2(dna, poly):
    return poly.count(TYR) > dna.getBeforeMutated().getPoly().count(TYR)

''' mutation 3 '''
mutationFrom3 = 'w'
mutationTo3 = 'z'
mutation3 = [
    (replace, DNA(1))
]
def condition3(dna, poly):
    return poly.count(GLN) > dna.getBeforeMutated().getPoly().count(GLN)