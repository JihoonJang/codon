from .DNA import DNA
from .Poly import Poly
from .Mutation import Mutation
from .Specification import *

class Wrapper:
    def __init__(self, dna, mutationList, condition, isPoly = False):
        self.iterator = None
        self.mutation = Mutation(dna, *mutationList[0])
        self.mutationList = mutationList
        self.dna = dna
        self.beforeMutated = []
        self.hashSet = set()
        for i in range(1, len(mutationList)):
            self.mutation = Mutation(self.mutation, *mutationList[i])
        self.condition = condition

    def __iter__(self):
        if self.iterator == None:
            self.mutation = iter(self.mutation)
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
                mut = next(self.mutation)
                mStack = mut.mutationStack
                mut.mutationStack = [mStack[i] for i in range(len(mStack) - len(self.mutationList), len(mStack))]
                if self.condition(mut, Poly(mut.getPoly())):
                    beforeMut = mut.getBeforeMutated()
                    h = hash(beforeMut)
                    if h not in self.hashSet:
                        self.hashSet.add(h)
                        self.beforeMutated.append(beforeMut)
                    return mut
        except StopIteration:
            raise StopIteration

    def getBeforeMutated(self):
        iter(self)
        return self.beforeMutated
    def beforeMutatedIsSameWith(self, iterator):
        iter(self)
        hashSet = set()
        for i in iterator:
            hashSet.add(hash(i))
        return hashSet == self.hashSet