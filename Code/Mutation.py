from .DNA import DNA
from .Function import allSame
from .Specification import stringToNucl

class Mutation:
    ''' dna : DNA type or Mutation type
        mutate : 돌연변이 종류 ('delete', 'insert', replace' 가능)
        mutateTo : 돌연변이 시 결실되는 or 삽입되는 or 치환 후의 염기
        mutateFrom : 치환 시 주형 가닥의 염기. 없으면 아무데서나 치환 가능
        ex) 1개의 퓨린 염기 결실 : Mutation(dna, 'delete', DNA('u'))
            2개의 연속된 동일한 염기 삽입 : Mutation(dna, 'insert', DNA(2, allSame))
            2개의 연속된 동일한 퓨린 염기가 각각 피리미딘 계열 염기로 치환 : Mutation(dna, 'replace', DNA('yy'), DNA('uu', allSame))

        만약 돌연변이가 2개 이상이면 다음과 같이 사용
        ex) 1개의 염기 결실 후 2개의 동일한 연속된 염기가 삽입
            m = Mutation(dna, 'delete', DNA(1))
            m = Mutation(m, 'insert', DNA(2, allSame))
            for dna in m:
                # do something
                # dna 중복 없음      
    '''
    def __init__(self, dna, mutate, mutateTo, mutateFrom = None):
        if 'DNA' in str(type(dna)):
            dna.mutationStack = []

        self.DNAList = dna

        self.mut = mutate
        self.mutToDNA = mutateTo
        
        if self.mut == 'delete' and self.mutToDNA.str == [stringToNucl['_']] * len(self.mutToDNA):
            self.deleteCount = len(self.mutToDNA)
        else:
            self.deleteCount = None

        if mutateFrom != None:
            self.mutFromDNA = mutateFrom
        else:
            self.mutFromDNA = None

    def __iter__(self):
        self.idx = 0
        self.DNAIterator = iter(self.DNAList)
        
        # StopIteration이 뜨면 안 됨
        try:
            self.curDNA = next(self.DNAIterator)
        except StopIteration:
            # 임시 방편
            self.curDNA = DNA('A')
        self.mutToDNAIterator = iter(self.mutToDNA)
        self.visit = set()
        
        return self

    def __next__(self):
        while True:
            try:
                mutDNA = next(self.mutToDNAIterator)
            except StopIteration:
                self.mutToDNAIterator = iter(self.mutToDNA)
                mutDNA = next(self.mutToDNAIterator)
                self.idx += 1
                if (self.mut == 'insert' and self.idx > len(self.curDNA)) or (self.mut != 'insert' and self.idx > len(self.curDNA) - len(mutDNA)):
                    try:
                        self.curDNA = next(self.DNAIterator)
                    except StopIteration:
                        raise StopIteration
                    self.idx = 0
            curDNA = self.curDNA.copy()
            if self.mut == 'insert':
                curDNA.insert(self.idx, mutDNA)
            elif self.mut == 'delete':
                if self.deleteCount != None:
                    curDNA.delete(self.idx, count = self.deleteCount)
                elif not curDNA.delete(self.idx, mutDNA):
                    continue
            elif self.mut == 'replace':
                if self.mutFromDNA == None:
                    if not curDNA.replace(self.idx, mutDNA):
                        continue
                else:
                    for fromDNA in self.mutFromDNA:
                        if curDNA.replace(self.idx, mutDNA, fromDNA):
                           break
                    else:
                        continue 
            else:
                print("mutation not exist (only delete, insert, replace possible)")
                assert False
            if hash(curDNA) not in self.visit:
                self.visit.add(hash(curDNA))
                return curDNA