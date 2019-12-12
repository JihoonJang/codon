from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *

''' 200915 예시 '''

standardDNA = DNA('CTATGCGGAGGATGGAAAGGAAGCTCTAGCTAG')


''' mutation 1 (w -> x)
    mutationFrom1 : w
    mutationTo1 : x
    dna[mutationFrom1] : mutationFrom1(w)의 염기 서열 (비주형 기준)
    mutation1 : 주형에서 CC 결실, 염기 하나 삽입
    condition1 : 6종류의 펩타이드, 3번째 아스파트산, 5번째 아르지닌
'''
mutationFrom1 = 'w'
mutationTo1 = 'x'
mutation1 = [
    (delete, DNA('CC')), 
    (insert, DNA(1))
]
def condition1(dna, poly):
    return poly.peptideKind() == 6 and poly.NthPeptide(3) == ASP and poly.NthPeptide(5) == ARG



''' mutation 2 (x -> y)
    mutation2 : 주형에서 T 결실, 염기 하나 삽입
    condition2 : 9종류의 펩타이드, 아스파트산과 히스티딘 가짐
'''
mutationFrom2 = 'x'
mutationTo2 = 'y'
mutation2 = [
    (delete, DNA('T')), 
    (insert, DNA(1))
]
def condition2(dna, poly):
    return poly.peptideKind() == 9 and poly.havePeptide(ASP) and poly.havePeptide(HIS)



''' mutation 3 (y -> z)
    mutation3 : 주형에서 2개 결실
    condition3 : 개수가 2개인 펩타이드의 종류가 2개 이상 (서로 다른 아미노산 ⓐ와 ⓑ를 각각 2개씩 가짐)
'''
mutationFrom3 = 'y'
mutationTo3 = 'z'
mutation3 = [
    (delete, DNA(2))
]
def condition3(dna, poly):
    return poly.peptideKind(2) >= 2