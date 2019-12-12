from Code.DNA import DNA
from Code.Function import *
from Code.Mutation import Mutation
from Code.Wrapper import Wrapper
from Condition import *


dnaStorage[mutationFrom1] = nonTemplateStrand

while True:
    try:
        dnaStorage[mutationTo1] = Wrapper(dnaStorage[mutationFrom1], mutation1, condition1)
    except:
        pass
    try:
        dnaStorage[mutationTo2] = Wrapper(dnaStorage[mutationFrom2], mutation2, condition2)
    except:
        pass
    try:
        dnaStorage[mutationTo3] = Wrapper(dnaStorage[mutationFrom3], mutation3, condition3)
    except:
        pass
    break




print(mutationFrom1, ":")
print("-" * 100)
for d in dnaStorage[mutationFrom1]:
    print(d)

''' mutation 1 '''

print(mutationTo1, ":")
print("-" * 100)
for d in dnaStorage[mutationTo1]:
    print(d)
    

''' mutation 2 '''
try:
    print(mutationTo2, ":")
    print("-" * 100)
    for d in dnaStorage[mutationTo2]:
        print(d)

except:
    pass

''' mutation 3 '''
try:
    print(mutationTo3, ":")
    print("-" * 100)
    for d in dnaStorage[mutationTo3]:
        print(d)

except:
    pass

print('end')