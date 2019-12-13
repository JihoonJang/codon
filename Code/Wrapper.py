from .DNA import DNA
from .Poly import Poly
from .Mutation import Mutation
from .Specification import *

class Wrapper:
    def __init__(self, dna, mutationList, condition):
        try:
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
            self.mutationToNormal = {}

            self.condition = condition
        except Exception as e:
            print(e)

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
                    if h not in self.mutationToNormal:
                        self.mutationToNormal[h] = beforeMut
                    else:
                        self.mutationToNormal[h] |= beforeMut
                    return mut
        except StopIteration:
            raise StopIteration

    def getBeforeMutated(self):
        iter(self)
        return list(self.mutationToNormal.values())

    def beforeMutatedIsDifferentWith(self, iterator):
        iter(self)
        hashSet = set()
        for i in iterator:
            hashSet.add(i.hash256())
        return not set(d.hash256() for d in self.getBeforeMutated()).issubset(hashSet)