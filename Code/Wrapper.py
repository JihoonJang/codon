from .DNA import DNA
from .Poly import Poly
from .Mutation import Mutation
from .Specification import *

class Wrapper:
    def __init__(self, dna, mutationList, condition):
        self.iterator = None
        
        for i in range(len(mutationList)):
            m = mutationList[i]
            new_m = (m[0], )
            for j in range(1, len(m)):
                if 'DNA' in str(type(m[j])):
                    new_m = (*new_m, m[j].complementReverse())
                elif j == 2:
                    new_m = (*new_m, None, m[j])
            mutationList[i] = new_m
            
                
        self.mutationList = mutationList
        self.dnaIter = dna
        self.beforeMutated = []
        self.filled = {}
        self.condition = condition
        # beforeFilledVisit : Mutation에서 중복 체크 없앤 후 추가 -> 아미노산 서열만 주어진 문제에서 실행 시간 단축 위해서
        self.beforeFilledVisit = {}
        self.visit = set()
    def __len__(self):
        iter(self)
        return len(self.iterator)
    def __iter__(self):
        if self.iterator == None:
            self.dnaIter = iter(self.dnaIter)
            try:
                self.dna = next(self.dnaIter)
            except StopIteration:
                return iter([])
            self.mutationIter = Mutation([self.dna], *self.mutationList[0])
            for i in range(1, len(self.mutationList)):
                self.mutationIter = Mutation(self.mutationIter, *self.mutationList[i])
            self.mutationIter = iter(self.mutationIter)
            self.iterator = []
            while True:
                try:
                    self.iterator.append(next(self))
                except StopIteration:
                    break
        iteratorHash = [d.hash256() for d in self.iterator]
        iterator = self.iterator
        self.iterator = []
        for i in range(len(iterator)):
            cnt = 0
            for h in iteratorHash:
                if h == h | iterator[i].hash256():
                    cnt += 1
            if cnt == 1:
                self.iterator.append(iterator[i])
                
        return iter(self.iterator)    
    
    def __next__(self):
        try:
            while True:
                # get mutation 
                try:
                    mut = next(self.mutationIter)
                except StopIteration:
                    self.dna = next(self.dnaIter)
                    self.mutationIter = Mutation([self.dna], *self.mutationList[0])
                    for i in range(1, len(self.mutationList)):
                        self.mutationIter = Mutation(self.mutationIter, *self.mutationList[i])
                    self.mutationIter = iter(self.mutationIter)
                    mut = next(self.mutationIter)   

                # mutation 초기화
                mStack = mut.mutationStack
                mut.mutationStack = [mStack[i] for i in range(len(mStack) - len(self.mutationList), len(mStack))]
                mut.poly = []

                # 염기서열 중복 체크 (테스트 필요)
                h1 = mut.hash256()
                if h1 in self.beforeFilledVisit:
                    self.beforeFilledVisit[h1].appendMutation(mut.mutationStack)
                    continue
                else:
                    mut.appendMutation(mut.mutationStack)
                    self.beforeFilledVisit[h1] = mut

                # 조건 만족 체크
                if self.condition(mut, Poly(mut.getPoly())):
                    if not mut.makeBlankFill():
                        continue
                    beforeMut = self.dna & mut.getBeforeMutated()
                    if beforeMut == None:
                        continue
                    beforeMut.poly = self.dna.poly
                    if not beforeMut.makeBlankFill():
                        continue
                    
                    h = self.dna.hash256()
                    if h not in self.filled:
                        self.filled[h] = beforeMut
                    else:
                        self.filled[h] |= beforeMut

                    h2 = mut.hash256()
                    if h2 not in self.visit:
                        self.visit.add(h2)
                        return mut
        except StopIteration:
            raise StopIteration

    def getBeforeMutated(self):
        iter(self)
        return list(self.filled.values())

    def beforeMutatedIsFilledThan(self, iterator):
        iter(self)
        hashSet = set()
        for i in iterator:
            hashSet.add(i.hash256())
        return not set(d.hash256() for d in self.getBeforeMutated()).issubset(hashSet)
    
    def __hash__(self):
        h = 0
        for i in self.iterator:
            h += i.hash256()
        return h