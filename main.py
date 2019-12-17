from Code.DNA import DNA
from Code.Function import *
from Code.Mutation import Mutation
from Code.Wrapper import Wrapper
from Condition import *
import time

startTime = time.time()
dnaStorage[mutationFrom1] = nonTemplateStrand

while True:
    try:
        dnaStorage[mutationTo1] = Wrapper(dnaStorage[mutationFrom1], mutation1, condition1)
        print(mutationFrom1, '->', mutationTo1, 'Calculating..')
        if dnaStorage[mutationTo1].beforeMutatedIsFilledThan(dnaStorage[mutationFrom1]):
            dnaStorage[mutationFrom1] = dnaStorage[mutationTo1].getBeforeMutated()
            continue
        dnaStorage[mutationFrom1] = dnaStorage[mutationTo1].getBeforeMutated()
    except NameError:
        pass
    try:
        dnaStorage[mutationTo2] = Wrapper(dnaStorage[mutationFrom2], mutation2, condition2)
        print(mutationFrom2, '->', mutationTo2, 'Calculating..')
        if dnaStorage[mutationTo2].beforeMutatedIsFilledThan(dnaStorage[mutationFrom2]):
            dnaStorage[mutationFrom2] = dnaStorage[mutationTo2].getBeforeMutated()
            continue
        if mutationFrom1 == mutationFrom2 and len(dnaStorage[mutationFrom2]) != len(dnaStorage[mutationTo2].getBeforeMutated()):
            dnaStorage[mutationFrom2] = dnaStorage[mutationTo2].getBeforeMutated()
            continue
        dnaStorage[mutationFrom2] = dnaStorage[mutationTo2].getBeforeMutated()
    except NameError:
        pass
    try:
        dnaStorage[mutationTo3] = Wrapper(dnaStorage[mutationFrom3], mutation3, condition3)
        print(mutationFrom3, '->', mutationTo3, 'Calculating..')
        if dnaStorage[mutationTo3].beforeMutatedIsFilledThan(dnaStorage[mutationFrom3]):
            dnaStorage[mutationFrom3] = dnaStorage[mutationTo3].getBeforeMutated()
            continue
        if mutationFrom2 == mutationFrom3 and len(dnaStorage[mutationFrom3]) != len(dnaStorage[mutationTo3].getBeforeMutated()):
            dnaStorage[mutationFrom3] = dnaStorage[mutationTo3].getBeforeMutated()
            continue
        dnaStorage[mutationFrom3] = dnaStorage[mutationTo3].getBeforeMutated()
    except NameError:
        pass
    break



print()
print(mutationFrom1, ":")
print("-" * 100)
for d in dnaStorage[mutationFrom1]:
    print(d)

''' mutation 1 '''
try:
    print(mutationTo1, ":")
    print("-" * 100)
    for d in dnaStorage[mutationTo1]:
        print(d)
except NameError:
    pass
    

''' mutation 2 '''
try:
    print(mutationTo2, ":")
    print("-" * 100)
    for d in dnaStorage[mutationTo2]:
        print(d)
except NameError:
    pass

''' mutation 3 '''
try:
    print(mutationTo3, ":")
    print("-" * 100)
    for d in dnaStorage[mutationTo3]:
        print(d)
except NameError:
    pass

print('End, Execution time :', time.time() - startTime, 'Sec')