from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *

''' 여기서 정상 주형 가닥 염기 서열, 돌연변이, 
    돌연변이 후 조건 수정 후 run 누르면 됨
    (main.py는 안 건드려도 됨!)
'''

''' 200915 예시 '''

''' w -> x 
    w : w의 염기 서열 (비주형 기준)
    w_mutation : 주형에서 CC 결실, 염기 하나 삽입
    w_cond : 6종류의 펩타이드, 3번째 아스파트산, 5번째 아르지닌
'''
w = DNA('CTATGCGGAGGATGGAAAGGAAGCTCTAGCTAG')
# w.makeCompRev()
w_mutation = [
    (delete, DNA('CC')), 
    (insert, DNA(1))
]
def w_cond(dna, poly):
    return poly.peptideKind() == 6 and poly.NthPeptide(3) == ASP and poly.NthPeptide(5) == ARG



''' x -> y 
    x_mutation : 주형에서 T 결실, 염기 하나 삽입
    x_cond : 9종류의 펩타이드, 아스파트산과 히스티딘 가짐
'''
x_mutation = [
    (delete, DNA('T')), 
    (insert, DNA(1))
]
def x_cond(dna, poly):
    return poly.peptideKind() == 9 and poly.havePeptide(ASP) and poly.havePeptide(HIS)



''' y -> z 
    y_mutation : 주형에서 2개 결실
    y_cond : 개수가 2개인 펩타이드의 종류가 2개 이상 (서로 다른 아미노산 ⓐ와 ⓑ를 각각 2개씩 가짐)
'''
y_mutation = [
    (delete, DNA(2))
]
def y_cond(dna, poly):
    return poly.peptideKind(2) >= 2