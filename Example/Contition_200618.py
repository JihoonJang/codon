from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *

''' 200618 예시 '''

standardDNA = DNA('메싸이오닌－발린－라이신－' + 가 + '－트레오닌－' + 나 + '－아이소류신－류신－글라이신')


''' mutation 1 '''
mutationFrom1 = 'x'
mutationTo1 = 'y'
가 = '아르지닌'
나 = '세린'
mutation1 = [
    (delete, DNA(1)),
    (insert, DNA(1))
]
def condition1(dna, poly):
    return dna.sequence('메싸이오닌－발린－세린－발린－히스티딘－글루타민－타이로신－발린－글라이신')


''' mutation 2 '''
mutationFrom2 = 'x'
mutationTo2 = 'z'
mutation2 = [
    (delete, DNA(2, allSame)), 
    (insert, DNA(2, allSame))
]
def condition2(dna, poly):
    return dna.mutationStack[0][1] & dna.mutationStack[1][1] != None and dna.polyLength(7) and poly.NthPeptide(4) | TYR > 0