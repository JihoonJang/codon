from .DNA import DNA
from .Function import allSame
from .Specification import stringToNucl
from .Position import Position

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
    '''
    def __init__(self, dna, mutate, mutateTo, mutateFrom = None, position = None):
        self.iterator = None

        self.DNAList = dna

        self.mut = mutate
        self.mutToDNA = mutateTo

        self.idxIterator = position
        
        if self.mut == 'delete' and self.mutToDNA.str == [stringToNucl['_']] * len(self.mutToDNA) and self.mutToDNA.condition == None:
            self.mutToDNA = [self.mutToDNA]

        self.mutFromDNA = mutateFrom

    def __iter__(self):
        if self.iterator == None:
            self.DNAIterator = iter(self.DNAList)
            
            # StopIteration이 뜨면 안 됨
            try:
                self.curDNA = next(self.DNAIterator)
            except StopIteration:
                return iter([])
            self.mutToDNAIterator = iter(self.mutToDNA)
            self.iterator = []

            if self.idxIterator == None:
                if self.mut == 'insert':
                    self.idxIterator = range(len(self.curDNA) + 1)
                elif self.mut == 'delete' or self.mut == 'replace':
                    self.idxIterator = range(len(self.curDNA) - len(self.mutToDNA) + 1)
                else:
                    print('mutate is not defined (not delete, insert, replace)')
                    assert False
            self.idxIteratorCur = iter(self.idxIterator)
            try:
                self.idx = next(self.idxIteratorCur)
            except StopIteration:
                print('mutate position error!')
                return iter([])
            while True:
                try:
                    self.iterator.append(next(self))
                except StopIteration:
                    break
        return iter(self.iterator)

    def __next__(self):
        while True:
            try:
                mutDNA = next(self.mutToDNAIterator)
            except StopIteration:
                self.mutToDNAIterator = iter(self.mutToDNA)
                mutDNA = next(self.mutToDNAIterator)
                try:
                    self.idx = next(self.idxIteratorCur)
                except StopIteration:
                    try:
                        self.curDNA = next(self.DNAIterator)
                    except StopIteration:
                        raise StopIteration
                    
                    self.idxIteratorCur = iter(self.idxIterator)
                    self.idx = next(self.idxIteratorCur)
            curDNA = self.curDNA.copy()
            if self.mut == 'insert':
                curDNA.insert(self.idx, mutDNA)
            elif self.mut == 'delete':
                if not curDNA.delete(self.idx, mutDNA):
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
            return curDNA