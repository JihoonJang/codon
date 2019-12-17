from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 190919 예시 '''

nonTemplateStrand = DNA('TTATGTTAGCTACCTTCCATCGTACGCATTAG'), DNA('TTATGTTAGCTACCTTCCATCGTACGCATTAG').complementReverse()


''' mutation 1 '''
mutationFrom1 = 'w'
mutationTo1 = 'x'
mutation1 = [
    (insert, DNA('u'))
]
def condition1(dna, poly):
    return (poly.length() == 4 or poly.length() == 9) and not poly.havePeptide(LEU) and not poly.havePeptide(SER)


''' mutation 2 '''
mutationFrom2 = 'w'
mutationTo2 = 'y'
mutation2 = [
    (insert, DNA('u'))
]
def condition2(dna, poly):
    return (poly.length() == 4 or poly.length() == 9) and not poly.havePeptide(LEU) and poly.count(SER) == 1 and poly.count(TYR) == 1

''' mutation 3 '''
mutationFrom3 = 'w'
mutationTo3 = 'z'
mutation3 = [
    (insert, DNA('u'))
]
def condition3(dna, poly):
    return (poly.length() == 4 or poly.length() == 9) and not poly.havePeptide(LEU) and poly.count(SER) == 2