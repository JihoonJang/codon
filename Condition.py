from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *

''' 여기서 정상 주형 가닥 염기 서열, 돌연변이, 
    돌연변이 후 조건 수정 후 run 누르면 됨
    (main.py는 안 건드려도 됨!)
'''

''' 200618 예시 '''

''' w -> x 
    w : w의 염기 서열 (비주형 기준)
    w_mutation : 주형에서 CC 결실, 염기 하나 삽입
    w_cond : 6종류의 펩타이드, 3번째 아스파트산, 5번째 아르지닌
'''
가 = '세린'
나 = '아르지닌'
x = DNA('메싸이오닌－발린－라이신－' + 가 + '－트레오닌－' + 나 + '－아이소류신－류신－글라이신')
# w.makeCompRev()
x_mutation = [
    (delete, DNA(1)), 
    (insert, DNA(1))
]
def x_cond(dna, poly):
    return dna.sequence('메싸이오닌－발린－세린－발린－히스티딘－글루타민－타이로신－발린－글라이신')



''' x -> y 
    x_mutation : 주형에서 T 결실, 염기 하나 삽입
    x_cond : 9종류의 펩타이드, 아스파트산과 히스티딘 가짐
'''
y_mutation = [
    (delete, DNA(2, allSame)), 
    (insert, DNA(2, allSame))
]
def y_cond(dna, poly):
    return dna.mutationStack[0][1].str == dna.mutationStack[1][1].str and poly.polyLength() == 6 and poly.NthPeptide() == TYR