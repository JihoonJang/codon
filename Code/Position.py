from .Specification import *
from .DNA import DNA

class Position:
    def __init__(self, dna, start, end = None, standard = 'DNA'):
        if end == None:
            end = start
        if standard != 'DNA':
            s = dna.getStartIdx()
            if s == -1:
                start = end = -1
            else:
                start = s + (start - 1) * 3
                end = s + end * 3 - 1
                
    def __iter__(self):
        return iter(range(start, end + 1))