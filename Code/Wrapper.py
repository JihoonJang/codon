from .DNA import DNA
from .Poly import Poly
from .Mutation import Mutation

class Wrapper:
    def __init__(self, dna, mutationList, condition):
        self.mutation = Mutation(dna, *mutationList[0])
        for i in range(1, len(mutationList)):
            self.mutation = Mutation(self.mutation, *mutationList[i])
        self.condition = condition
    
    def __iter__(self):
        self.mutation = iter(self.mutation)
        return self
    
    def __next__(self):
        try:
            while True:
                mut = next(self.mutation)
                if self.condition(mut, Poly(mut.getPoly())):
                    return mut
        except StopIteration:
            raise StopIteration
