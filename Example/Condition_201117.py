from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *
from itertools import permutations

''' 201117 예시 '''

standardDNA = DNA('GACTCACAAGCCATTGAACCAACTCGTTGCCATGC').complementReverse()

''' mutation 1 (y -> x)
    mutationFrom1 : y
    mutationTo1 : x
    dna[mutationFrom1] : mutationFrom1(w)의 염기 서열 (비주형 기준)
    mutation1 : 주형에서 연속된 2개 염기 삽입, 연속된 2개의 염기 결실
    condition1 : MET-ALA-Seq[0]-Seq[1]-Seq[2]의 아미노산 서열을 가짐 (Seq[0], Seq[1], Seq[2] : perm 참조)
'''
mutationFrom1 = 'y'
mutationTo1 = 'x'

mutation1 = [
    (delete, DNA(2)), 
    (insert, DNA(2))
]
def condition1(dna, poly):
    perm = permutations(['류신-발린', '발린-글루타민-트립토판', '라이신-류신'])
    for seq in perm:
        if poly.haveSequence('메싸이오닌-알라닌-' + seq[0] + '-' + seq[1] + '-' + seq[2]):
            return True
    return False



''' mutation 2 (x -> z)
    mutation2 : 주형에서 C 결실
    condition2 : 6종류 아미노산, 4번째 트립토판
'''
mutationFrom2 = 'x'
mutationTo2 = 'z'
mutation2 = [
    (delete, DNA('C'))
]
def condition2(dna, poly):
    return poly.peptideKind() == 6 and poly.NthPeptide(4) == TRP