from .Specification import *

class DNA:
    ''' Attributes
            str : DNA 염기 서열 저장 (비주형 기준!!)
            poly : 염기 서열 대신 polypeptide 서열만 주어졌을 경우 사용됨 (ex : 160618)
            mutationStack : 돌연변이 내역을 stack에 저장
    '''
    def __init__(self, arg = None, condition = None):
        self.poly = []
        self.str = []
        self.mutationStack = []
        self.condition = condition
        if type(arg) == type(''):
            if arg[0] in stringToNucl:
                self.str = []
                for n in arg:
                    self.str.append(stringToNucl[n])
            else:
                self.poly = stringToPoly(arg)
                self.poly.append(END)
                for p in self.poly:
                    try:
                        for n in peptideToTriplet[p]:
                            self.str.append(n)
                    except KeyError:
                        for i in range(3):
                            self.str.append(stringToNucl['_'])
                self.poly.pop()
        elif type(arg) == type(0):
            self.str = [stringToNucl['_']] * arg

    def printDNA(self):
        print('5－', end='')
        for n in self.str:
            assert n != 0
            print(nuclToString[n], end='')
        print('－3')

    def printPoly(self):
        poly = self.getPoly()
        if len(poly) == 0:
            print("합성 안 됨")
            return
        for p in poly[0:len(self.poly) - 1]:
            try:
                print(peptideToString[p], end='－')
            except:
                if p & END == 0:
                    print('___', end='－')
                else:
                    print('(종결)', end='－')
        try:
            print(peptideToString[poly[-1]])
        except:
            if poly[-1] & END == 0:
                print('___')
            else:
                print('(종결)')
    
    def printMutation(self):
        for m in self.mutationStack:
            print(m[0], end=', 위치 : ')
            if m[0] == 'replace':
                print(m[3], end = ', 염기 서열 : ')
            else:
                print(m[2], end = ', 염기 서열 : ')
            m[1].printDNA()

    def copy(self, polyCopy = True):
        ''' return : copy of self '''
        d = DNA()
        d.str = self.str.copy()
        if polyCopy:
            d.poly = self.poly.copy()
        d.mutationStack = self.mutationStack.copy()
        d.condition = self.condition
        return d
    
    def makeCompRev(self):
        compRev = []
        self.str.reverse()
        for n in self.str:
            newN = 0
            if n & A > 0:
                newN |= T
            if n & T > 0:
                newN |= A
            if n & G > 0:
                newN |= C
            if n & C > 0:
                newN |= G
            compRev.append(newN)
        self.str = compRev
        
    def __len__(self):
        ''' ex) len(DNA('AATTTGC')) == 7 '''
        return len(self.str)

    # 다른 attribute도 복제해야 하는지?
    def __getitem__(self, i):
        ''' i : index (1, 1:2, 1:5:2 등 list의 slicing과 같이 사용 가능)
            return : subDNA (type : DNA), copy를 반환
            ex) DNA('AATTTGC')[2:4] == DNA('TTT')
        '''
        try:
            d = DNA()
            if type(i) == type(0):
                d.str = [self.str[i]]
            else:
                d.str = self.str[i]
            return d
        except Exception as e:
            print('__getitem__ error : ', e)
    
    def __hash__(self):
        return hash(tuple(self.str))

    # DNA1 | DNA2 : 길이가 같지 않으면 None return, 같으면 각 원소에 대해 or 연산 수행 후 return
    def __or__(self, dna):
        try:
            if type(dna) == type(''):
                dna = DNA(dna)
            if len(self) != len(dna):
                return None
            ret = DNA(len(self))
            for i in range(len(self)):
                ret.str[i] = self.str[i] | dna.str[i]
            return ret
        except:
            assert('__or__ error')
    
    # DNA1 & DNA2 : 길이가 같지 않거나 and 연산 결과가 0이면 None return, 같으면 각 원소에 대해 and 연산 수행 후 return
    def __and__(self, dna):
        try:
            if type(dna) == type(''):
                dna = DNA(dna)
            if len(self) != len(dna):
                return None
            ret = DNA(len(self))
            for i in range(len(self)):
                ret.str[i] = self.str[i] & dna.str[i]
                if ret.str[i] == 0:
                    return None
            return ret
        except Exception as e:
            print('__and__ error : ', e)

    def __iter__(self):
        self.iterateIdx = []
        self.iterateNucl = []
        self.iterateDNA = self.copy(False)
        self.iterateEnd = False
        s = self.str

        if len(self.poly) == 0:
            for i in range(len(self) - 1, -1, -1):
                if s[i] == A or s[i] == T or s[i] == G or s[i] == C:
                    continue
                
                nucl = A
                while nucl != NUCL_END:
                    if s[i] & nucl > 0:
                        break
                    nucl <<= 1
                assert nucl != NUCL_END
                self.iterateIdx.append(i)
                self.iterateNucl.append(nucl)

            if len(self.iterateNucl):
                self.iterateNucl[0] = 0
            for i in range(len(self.iterateIdx)):
                self.iterateDNA.str[self.iterateIdx[i]] = self.iterateNucl[i]
        return self
        
    def __next__(self):
        def getNextNucl(i, curNucl):
            if curNucl == 0:
                nextNucl = A
            else:
                nextNucl = curNucl << 1
            while nextNucl != NUCL_END:
                if nextNucl & self.str[i] > 0:
                    break
                nextNucl <<= 1
            if nextNucl == NUCL_END:
                return None
            return nextNucl
        
        def getNext(self):
            self.iterateDNA = self.iterateDNA.copy(False)
            if len(self.iterateIdx) == 0 or len(self.poly) > 0:
                if self.iterateEnd:
                    raise StopIteration
                self.iterateEnd = True
                return self.iterateDNA
            for i in range(len(self.iterateIdx)):
                curIdx = self.iterateIdx[i]
                nextNucl = getNextNucl(curIdx, self.iterateNucl[i])
                if nextNucl != None:
                    self.iterateDNA.str[curIdx] = self.iterateNucl[i] = nextNucl
                    break
            else:
                raise StopIteration

            for j in range(i):
                curIdx = self.iterateIdx[j]
                self.iterateDNA.str[curIdx] = self.iterateNucl[j] = getNextNucl(curIdx, 0)
        
        getNext(self)
        while self.condition != None and self.condition(self.iterateDNA.str) == False:
            getNext(self)
        return self.iterateDNA 

    def insert(self, offset, strand):
        ''' offset : 시작 index
            strand : 삽입되는 DNA (type : DNA or str)
            return : None
            ex) DNA('AATTTGC').insert(2, 'GG') : DNA('AAGGTTTGC')가 됨
        '''
        try:
            strand.makeCompRev()
            if type(strand) == type(''):
                strand = DNA(strand)
            self.str[offset:offset] = strand.str
            self.mutationStack.append(('insert', strand, offset))
        except Exception as e:
            print("insertError : ", e)

    def delete(self, offset, strand = None, count = 0):
        ''' offset : 시작 index
            count : 결실시키는 연속된 염기의 개수
            strand : 지우고자 하는 DNA 염기 (type : DNA or str), 170618같은 문제에서만 사용
            return : False if strand != None and subDNA & strand == None else True
            ex) DNA('AATTTGC').delete(2, 3) == True, DNA('AAGC')가 됨
                DNA('AATTTGC').delete(2, 'TTT') == True, DNA('AAGC')가 됨
                DNA('AATTTGC').delete(2, DNA('333')) == False
        '''
        try:
            if strand == None:
                self.mutationStack.append(('delete', self[offset:offset+count], offset))
                self.str[offset:offset + count] = []
            else:
                if type(strand) == type(''):
                    strand = DNA(strand)
                strand.makeCompRev()
                d = self[offset:offset+len(strand)] & strand
                if d == None:
                    return False
                self.mutationStack.append(('delete', d, offset))
                self.str[offset:offset + len(strand)] = []
            return True
        except Exception as e:
            print("deleteError : ", e)

    def replace(self, offset, strandTo, strandFrom = None):
        ''' offset : 시작 index
            strandTo : 치환시키는 염기 서열
            strandFrom : 주형 가닥에서 치환시키고자 하는 연속된 염기 서열
            return : False if (strandFrom != None and subDNA & strandFrom == None) or strandTo == strandFrom else True
            ex) DNA('AATTTGC').replace(2, 'GGG') == True, DNA('AAGGGGC')가 됨
                DNA('AA222GC').replace(2, 'GGG', 'TTT') == True, DNA('AAGGGGC')가 됨
                DNA('AATTTGC').replace(2, 'TTT') == False
                DNA('AA222GC').replace(2, 'GGG', 'CCC') == False
        '''
        try:
            if type(strandTo) == type(''):
                strandTo = DNA(strandTo)
            strandTo.makeCompRev()
            if strandFrom == None:
                strandFrom = self[offset:offset+len(strandTo)]
                for i in range(len(strandFrom.str)):
                    if strandFrom.str[i] == strandTo.str[i]:
                        return False
                self.str[offset:offset+len(strandTo)] = strandTo.str
                self.mutationStack.append(('replace', strandTo, strandFrom, offset))
            else:
                if type(strandFrom) == type(''):
                    strandFrom = DNA(strandFrom)
                strandFrom.makeCompRev()
                d = self[offset:offset+len(strandFrom)] & strandFrom
                if d == None:
                    return False
                for i in range(len(strandFrom.str)):
                    if strandFrom.str[i] == strandTo.str[i]:
                        return False
                self.str[offset:offset+len(strandFrom)] = strandTo.str
                self.mutationStack.append(('replace', strandTo, d, offset))
            return True
        except Exception as e:
            print("replaceError : ", e)

    def getPoly(self):
        def _tripletToPeptide(tri, endFlag):
            try:
                return tripletToPeptide[tuple(tri.str)]
            except:
                endFlag[0] = False
                pep = 0
                for t in tri:
                    pep |= tripletToPeptide[tuple(t.str)]
                return pep
        try:
            endFlag = [True]
            if len(self.poly) > 0:
                return self.poly
            poly = []
            for i in range(len(self.str) - 2):
                if _tripletToPeptide(self[i : i + 3], endFlag) & MET > 0:
                    break
            else:
                return []
            
            for j in range(i, len(self.str) - 2, 3):
                poly.append(_tripletToPeptide(self[j : j + 3], endFlag))
                if poly[-1] == END > 0:
                    poly.pop()
                    break
            else:
                if endFlag[0]:
                    return []
            return poly
        except Exception as e:
            print("getPolyError : ", e)

    def getList(self, isSame = False):
        try:
            def getNext(base, curBase):
                originBase = base & curBase
                if curBase == 0:
                    curBase = 1
                    if base & curBase != 0:
                        return curBase
                while True:
                    curBase *= 2
                    if curBase >= NUCL_END:
                        return originBase
                    if base & curBase != 0:
                        return curBase
                return originBase
                
            result = []
            offset = self.str.copy()
            for i in range(len(offset)):
                offset[i] = getNext(offset[i], 0)
            result.append(offset)

            while True:
                new = self.str.copy()
                for i in range(len(new)):
                    new[i] = getNext(new[i], offset[i])
                if new == offset:
                    return result
                result.append(new)
                offset = new
        except:
            print('getListError')

    def sequence(self, pp):
        pp = DNA(pp)
        tmp = (pp & self)
        if tmp != None:
            self.str = tmp.str
            return True
        else:
            return False
    
    def polyLength(self, n):
        try:
            if self.getPoly()[n - 1] | END > 0:
                return True
            return False
        except:
            return False