from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 180617 예시 '''
nonTemplateStrand = DNA('CTATGCTGCATGGACGTTGCGACCGACCATAGGAT')


''' mutation 1 '''
mutationFrom1 = 'y'
mutationTo1 = 'x'
mutation1 = [
    (delete, DNA(1)),
    (insert, DNA(1))
]
def condition1(dna, poly):
    return poly.count(VAL) == 2 and poly.count(ALA) == 1 and poly.count(ASP) == 1 and poly.count(PRO) == 1 and poly.count(LEU) == 1 and poly.count(THR) == 1 and poly.count(HIS) == 1

''' mutation 2 '''
mutationFrom2 = 'x'
mutationTo2 = 'z'
mutation2 = [
    (delete, DNA(2, allSame)),
    (insert, DNA(2, allSame))
]
def condition2(dna, poly):
    return poly.haveSequence('메싸이오닌-류신-아스파라진-메싸이오닌-트레오닌-류신-아르지닌-프롤린')