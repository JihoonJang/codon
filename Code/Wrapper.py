from .DNA import DNA
from .Poly import Poly
from .Mutation import Mutation
from .Specification import *

class Wrapper:
    def __init__(self, dna = None, mutationList = None, condition = None, iterator = None):
        if iterator == None:
            self.iterator = None
            for i in range(len(mutationList)):
                m = mutationList[i]
                if len(m) == 2:
                    mutationList[i] = (m[0], m[1].complementReverse())
                elif len(m) == 3:
                    mutationList[i] = (m[0], m[1].complementReverse(), m[2].complementReverse())
            self.mutationList = mutationList
            self.dnaIter = dna
            self.beforeMutated = []
            self.filled = {}
            self.condition = condition
            self.visit = set()

        else:
            self.iterator = iterator

    def __iter__(self):
        if self.iterator == None:
            self.dnaIter = iter(self.dnaIter)
            self.dna = next(self.dnaIter)
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
                try:
                    mut = next(self.mutationIter)
                except StopIteration:
                    self.dna = next(self.dnaIter)
                    self.mutationIter = Mutation([self.dna], *self.mutationList[0])
                    for i in range(1, len(self.mutationList)):
                        self.mutationIter = Mutation(self.mutationIter, *self.mutationList[i])
                    self.mutationIter = iter(self.mutationIter)
                    mut = next(self.mutationIter)     
                mStack = mut.mutationStack
                mut.mutationStack = [mStack[i] for i in range(len(mStack) - len(self.mutationList), len(mStack))]
                mut.poly = []
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
                    h = mut.hash256()
                    if h not in self.visit:
                        self.visit.add(h)
                        return mut
        except StopIteration:
            raise StopIteration

    def getBeforeMutated(self):
        iter(self)
        return Wrapper(iterator = list(self.filled.values()))

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