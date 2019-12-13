from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 170618 예시 '''
nonTemplateStrand = DNA('메싸이오닌-메싸이오닌-아르지닌-트립토판-트레오닌-류신-글루타민-알라닌-아이소류신')


''' mutation 1 '''
mutationFrom1 = 'x'
mutationTo1 = 'x*'
mutation1 = [
    (delete, DNA(1)),
    (insert, DNA(1))
]
def condition1(dna, poly):
    return dna.sequenceIs('메싸이오닌-메싸이오닌-아르지닌-세린-아스파트산-발린-알라닌-트레오닌-아이소류신')

''' mutation 2 '''
mutationFrom2 = 'x'
mutationTo2 = 'x**'
mutation2 = [
    (delete, DNA(2, allSame)),
    (insert, DNA(2, allSame))
]
def condition2(dna, poly):
    return dna.sequenceIs('메싸이오닌-아이소류신-세린-아스파트산-글라이신-(가)-글루타민-알라닌-아이소류신')