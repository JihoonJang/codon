A, T, G, C, NUCL_END = [(1 << i) for i in range(5)]

PHE, LEU, SER, TYR, CYS, TRP, PRO,  \
HIS, GLN, ARG, IIE, MET, THR, ASN,  \
LYS, VAL, ALA, ASP, GLU, GLY, END,   \
PEPTIDE_END = [(1 << i) for i in range(22)]
PEPTIDE_ALL = PEPTIDE_END - 1

tripletToPeptide = {
    (T, T, T) : PHE, (T, T, C) : PHE, (T, T, A) : LEU, (T, T, G) : LEU,
    (T, C, T) : SER, (T, C, C) : SER, (T, C, A) : SER, (T, C, G) : SER,
    (T, A, T) : TYR, (T, A, C) : TYR, (T, A, A) : END, (T, A, G) : END,
    (T, G, T) : CYS, (T, G, C) : CYS, (T, G, A) : END, (T, G, G) : TRP,
    (C, T, T) : LEU, (C, T, C) : LEU, (C, T, A) : LEU, (C, T, G) : LEU,
    (C, C, T) : PRO, (C, C, C) : PRO, (C, C, A) : PRO, (C, C, G) : PRO,
    (C, A, T) : HIS, (C, A, C) : HIS, (C, A, A) : GLN, (C, A, G) : GLN,
    (C, G, T) : ARG, (C, G, C) : ARG, (C, G, A) : ARG, (C, G, G) : ARG,
    (A, T, T) : IIE, (A, T, C) : IIE, (A, T, A) : IIE, (A, T, G) : MET,
    (A, C, T) : THR, (A, C, C) : THR, (A, C, A) : THR, (A, C, G) : THR,
    (A, A, T) : ASN, (A, A, C) : ASN, (A, A, A) : LYS, (A, A, G) : LYS,
    (A, G, T) : SER, (A, G, C) : SER, (A, G, A) : ARG, (A, G, G) : ARG,
    (G, T, T) : VAL, (G, T, C) : VAL, (G, T, A) : VAL, (G, T, G) : VAL,
    (G, C, T) : ALA, (G, C, C) : ALA, (G, C, A) : ALA, (G, C, G) : ALA,
    (G, A, T) : ASP, (G, A, C) : ASP, (G, A, A) : GLU, (G, A, G) : GLU,
    (G, G, T) : GLY, (G, G, C) : GLY, (G, G, A) : GLY, (G, G, G) : GLY
}

peptideToString = {
    PHE : '페닐알라닌',
    LEU : '류신',
    SER : '세린',
    TYR : '타이로신',
    END : '종결',
    CYS : '시스테인',
    TRP : '트립토판',
    PRO : '프롤린',
    HIS : '히스티딘',
    GLN : '글루타민',
    ARG : '아르지닌',
    IIE : '아이소류신',
    MET : '메싸이오닌',
    THR : '트레오닌',
    ASN : '아스파트산',
    LYS : '라이신',
    VAL : '발린',
    ALA : '알라닌',
    ASP : '아스파트산',
    GLU : '글루탐산',
    GLY : '글라이신'
}

stringToPeptide = {v: k for k, v, in peptideToString.items()}


stringToNucl = {
    'A' : A,
    'T' : T,
    'G' : G,
    'C' : C,
    'u' : A | G,
    'y' : T | C,
    'k' : A | C,
    'i' : A | T | C,
    'j' : T | G,
    '2' : A | T,
    '3' : G | C,
    '_' : A | T | G | C
}

nuclToString = {v: k for k, v, in stringToNucl.items()}

peptideToTriplet = {}

delete = 'delete'
insert = 'insert'
replace = 'replace'

def stringToPoly(pp):
    poly = pp.replace('-', '－').split('－')
    for i in range(len(poly)):
        try:
            poly[i] = stringToPeptide[poly[i]]
        except KeyError:
            poly[i] = PEPTIDE_ALL
    return poly

dnaStorage = {}
polyFlag = False
seqChangeFlag = True

def printWithIndent(indent, *args):
    for s in args:
        s = str(s).split('\n')
        for i in s:
            if indent > 0:
                print(' ' * indent, i)
            else:
                print(i)