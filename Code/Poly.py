from .Specification import *

class Poly:
    ''' poly : list type 
        ex) self.poly = [MET, THR, SER, LEU, ASP, GLY, TRP]
    '''
    def __init__(self, poly):
        self.poly = poly

    def peptideKind(self, count = None):
        # 몇 종류의 펩타이드인지 return
        if count == None:
            return len(set(self.poly))
        else:
            dic = dict()
            for p in self.poly:
                if p in dic:
                    dic[p] += 1
                else:
                    dic[p] = 1
            ret = 0
            for key in dic:
                if dic[key] == count:
                    ret += 1
            return ret
    
    def length(self):
        # 폴리펩타이드의 길이 return
        return len(self.poly)
    
    def peptideCount(self, peptide):
        # peptide가 몇 개 있는지 return
        return self.poly.count(peptide)

    def NthPeptide(self, i):
        # n번째 펩타이드 return
        try:
            return self.poly[i - 1]
        except:
            return 0

    def havePeptide(self, peptide):
        # peptide를 가지면 True
        return peptide in self.poly

    def haveSequence(self, subSeq):
        # 특정 부분 서열을 가지면 True
        pp = stringToPoly(subSeq)
        for i in range(len(self.poly) - len(pp) + 1):
            if pp == self.poly[i:i+len(pp)]:
                return True
        return False