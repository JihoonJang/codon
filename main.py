from Code.DNA import DNA
from Code.Function import *
from Code.Mutation import Mutation
from Code.Wrapper import Wrapper
from Condition import *


try:
    print('w : ')
    w.printDNA()
    w.printPoly()
    print()
    x = Wrapper(w, w_mutation, w_cond)
except:
    x = x


print("x :")
for pp in x:
    pp.printDNA()
    pp.printPoly()
    print()

y = Wrapper(x, x_mutation, x_cond)
print("y : ")
for pp in y:
    pp.printDNA()
    pp.printPoly()
    print()

z = Wrapper(y, y_mutation, y_cond)
print("z : ")
for pp in z:
    pp.printDNA()
    pp.printPoly()
    print()