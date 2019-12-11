from Code.DNA import DNA
from Code.Function import *
from Code.Mutation import Mutation
from Code.Wrapper import Wrapper
from Condition import *



print(mutationFrom1, ":")
print("-" * 100)
for d in dna[mutationFrom1]:
    d.printDNA()
    d.printPoly()
    print()

''' mutation 1 '''
dna[mutationTo1] = Wrapper(dna[mutationFrom1], mutation1, condition1)

print(mutationTo1, ":")
print("-" * 100)
for d in dna[mutationTo1]:
    d.printDNA()
    d.printPoly()
    d.printMutation()
    print()
    

''' mutation 2 '''
try:
    dna[mutationTo2] = Wrapper(dna[mutationFrom2], mutation2, condition2)

    print(mutationTo2, ":")
    print("-" * 100)
    for d in dna[mutationTo2]:
        d.printDNA()
        d.printPoly()
        d.printMutation()
        print()

except:
    print(end='')


''' mutation 3 '''
try:
    dna[mutationTo3] = Wrapper(dna[mutationFrom3], mutation3, condition3)

    print(mutationTo3, ":")
    print("-" * 100)
    for d in dna[mutationTo3]:
        d.printDNA()
        d.printPoly()
        d.printMutation()
        print()

except:
    print(end='')

print('end')