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
        self.iterator = None
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

    def getDNAStr(self):
        s = '5－'
        for n in self.str:
            assert n != 0
            s += nuclToString[n]
        s += '－3'
        return s

    def getPolyStr(self):
        poly = self.getPoly()
        s = ''
        if len(poly) == 0:
            s = '합성 안 됨'
            return s
        for p in poly[0:len(self.poly) - 1]:
            try:
                s += peptideToString[p] + '－'
            except:
                if p & END == 0:
                    s += '___－'
                else:
                    s += '(종결)－'
        try:
            s += peptideToString[poly[-1]]
        except:
            if poly[-1] & END == 0:
                s += '___'
            else:
                s += '(종결)'
        return s

    
    def getMutationStr(self):
        s = ''
        for m in self.mutationStack:
            s += m[0] + ', 위치 : '
            if m[0] == 'replace':
                s += str(m[3]) + ', 염기 서열 : '
            else:
                s += str(m[2]) + ', 염기 서열 : '
            s += m[1].getDNAStr() + '\n'
        return s
    
    def __str__(self):
        s = self.getDNAStr() + '\n'
        s += self.getPolyStr() + '\n'
        s += self.getMutationStr() + '\n'
        return s

    def copy(self):
        ''' return : copy of self '''
        d = DNA()
        d.str = self.str.copy()
        d.poly = self.poly.copy()
        d.mutationStack = self.mutationStack.copy()
        d.condition = self.condition
        try:
            d.iterator = self.iterator.copy()
        except AttributeError:
            d.iterator = None
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
    
    def complementReverse(self):
        ret = self.copy()
        ret.makeCompRev()
        return ret

    def getBeforeMutated(self):
        d = self.copy()
        mStack = self.mutationStack.copy()
        while len(mStack) > 0:
            m = mStack.pop()
            if m[0] == 'delete':
                d.insert(m[2], m[1])
            elif m[0] == 'insert':
                d.delete(m[2], count = len(m[1]))
            else:
                d.replace(m[3], m[2])
        d.mutationStack = []
        return d
        
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

    def hash256(self):
        h = 0
        i = 0
        for n in self.str:
            h += n * (1 << i)
            i += 4
        return h

    # DNA1 | DNA2 : 길이가 같지 않으면 None return, 같으면 각 원소에 대해 or 연산 수행 후 return
    def __or__(self, dna):
        try:
            if type(dna) == type(''):
                dna = DNA(dna)
            if len(self) != len(dna):
                return None
            ret = self.copy()
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
            ret = self.copy()
            for i in range(len(self)):
                ret.str[i] = self.str[i] & dna.str[i]
                if ret.str[i] == 0:
                    return None
            return ret
        except Exception as e:
            print('__and__ error : ', e)

    def __iter__(self):
        if self.iterator == None:
            self.iterateDNA = self.copy()
            if len(self.poly) > 0:
                self.iterator = [self.iterateDNA]
                return iter(self.iterator)
            self.iterateIdx = []
            self.iterateNucl = []
            self.iterateEnd = False
            s = self.str

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

            if len(self.iterateIdx) == 0:
                self.iterator = [self.iterateDNA]
                return iter(self.iterator)
            
            self.iterator = []
            while True:
                try:
                    self.iterator.append(next(self))
                except StopIteration:
                    break
        return iter(self.iterator)
        
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
            self.iterateDNA = self.iterateDNA.copy()
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
            if type(strand) == type(''):
                strand = DNA(strand)
            self.str[offset:offset] = strand.str
            self.mutationStack.append(('insert', strand.copy(), offset))
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
            except KeyError:
                endFlag[0] = False
                pep = 0
                for t in tri:
                    pep |= tripletToPeptide[tuple(t.str)]
                return pep
        try:
            endFlag = [True]
            poly = []
            for i in range(len(self.str) - 2):
                if _tripletToPeptide(self[i : i + 3], endFlag) & MET > 0:
                    break
            else:
                return []
            
            for j in range(i, len(self.str) - 2, 3):
                poly.append(_tripletToPeptide(self[j : j + 3], endFlag))
                if poly[-1] & END > 0:
                    poly.pop()
                    break
            else:
                if endFlag[0]:
                    return []
            return poly
        except Exception as e:
            print("getPolyError : ", e)

    def makeBlankFill(self):
        try:
            if len(self.poly) == 0:
                return True
            for i in range(len(self.poly)):
                if not self.NthPeptideIs(i + 1, self.poly[i]):
                    return False
            try:
                return self.NthPeptideIs(len(self.poly) + 1, END)
            except:
                return False
        except Exception as e:
            print('makeBlankFill error :', e)

    def sequence(self, pp):
        try:
            pp = DNA(pp)
            tmp = pp & self
            if tmp != None:
                self.str = tmp.str
                self.poly = pp.poly
                return True
            else:
                return False
        except Exception as e:
            print('DNA sequence error :', e)
    
    def polyLength(self, n):
        try:
            if self.getPoly()[n - 1] | END > 0:
                return True
            return False
        except:
            return False

    def NthPeptideIs(self, i, p):
        i = (i - 1) * 3
        if p == SER:
            if self.str[i + 2] & (T | C) == 0:  # __u인 경우 -> TCu만 가능 (AGu는 불가)
                self.str[i] = T
                self.str[i + 1] = C
            elif self.str[i] == T:
                self.str[i + 1] = C
            elif self.str[i] == A:
                self.str[i + 1] = G
                self.str[i + 2] &= T | C
            elif self.str[i + 1] == C:
                self.str[i] = T
            elif self.str[i + 1] == G:
                self.str[i] = A
                self.str[i + 2] &= T | C
        elif p == LEU:
            if self.str[i] == T:
                self.str[i + 2] &= A | G
            elif self.str[i + 2] & (A | G) == 0:
                self.str[i] = C
        elif p == ARG:
            if self.str[i] == A:
                self.str[i + 2] &= A | G
            elif self.str[i + 2] & (A | G) == 0:
                self.str[i] = C
        elif p == END:
            if self.str[i + 1] == G and self.str[i + 2] == G:
                return False
        if p in peptideToTriplet:
            for j in range(3):
                self.str[i + j] &= peptideToTriplet[p][j]
                if self.str[i + j] == 0:
                    return False
        else:
            # peptide가 주어지지 않은 경우
            pass
        return True