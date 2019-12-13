from .Specification import *

for tri in tripletToPeptide:
    pep = tripletToPeptide[tri]
    if pep in peptideToTriplet:
        peptideToTriplet[pep] = tuple([tri[i] | peptideToTriplet[pep][i] for i in range(3)])
    else:
        peptideToTriplet[pep] = tri