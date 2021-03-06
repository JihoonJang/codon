# 사용법
----
##### 수정할 부분
Condition.py 파일을 문제 조건에 맞게 변경한 후 run을 누르면 됩니다. 다른 파일은 수정하지 않아도 됩니다.


```python
''' mutation (w -> x) '''
mutationFrom = 'w' 
mutationTo = 'x'   

dna[mutationFrom] = DNA('ATGGTACGTACGATACGACGACGTTGACA')
mutation = [
    (delete, DNA(2, allSame)),                  # 2개의 동일한 (allSame) 염기 결실
    (replace, DNA('yy'), DNA('uu', allSame))    # 2개의 연속된 동일한 퓨린 염기('uu')가 각각 피리미딘 염기('yy')로 치환
]
def condition(dna, poly):
    # 6개의 서로 다른 아미노산으로 구성 (= 6개의 아미노산을 가지고, 펩타이드의 종류가 6개)
    return poly.length() == 6 and poly.peptideKind() == 6
```

위 코드에서 `mutateFrom`, `mutateTo`, `dna[mutationFrom]`, `mutation`, `condition`을 문제 조건에 맞게 변형해서 사용하시면 됩니다.

##### 예시 (Contition_200915.py)
```python
from Code.DNA import DNA
from Code.Poly import Poly
from Code.Function import allSame
from Code.Specification import *

''' 200915 예시 '''

''' mutation 1 (w -> x)
    mutationFrom1 : w
    mutationTo1 : x
    dna[mutationFrom1] : w의 염기 서열 (비주형 기준)
    mutation1 : 주형에서 CC 결실, 염기 하나 삽입
    condition1 : 6종류의 펩타이드, 3번째 아스파트산, 5번째 아르지닌
'''
mutationFrom1 = 'w'
mutationTo1 = 'x'
dna[mutationFrom1] = DNA('CTATGCGGAGGATGGAAAGGAAGCTCTAGCTAG')
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
```
