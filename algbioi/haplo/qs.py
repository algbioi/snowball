"""
    Copyright (C) 2015  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Module that handle quality scores of consensus sequences. Module version 1.2
"""

import numpy as np

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

PROB_MATRIX_ROWS = 5  # rows correspond to: A, C, T, G, coverage


def nuclToInt(nucl):
    """
        Map: {A/a, C/c, T/t/U/u, G/g} -> {0, 1, 2, 3}
        Make 'xor 0b10" on the return value to get the reverse complement nucleotides.

        @raise TypeError: for non-DNA characters (i.e. not included in the mapping)
        @type nucl: str
        @rtype: int
    """
    i = ord(nucl)
    # A/a/C/c or G/g or T/t/U/u
    if (i | 0b1100010 == 0b1100011) or (i | 0b1100000 == 0b1100111) or (i | 0b1100001 == 0b1110101):
        return (i >> 1) & 0b11
    else:
        raise TypeError('intToNucl: argument "%s" not allowed' % nucl)


def isDNAChar(ch):
    """
        @type ch: str
        @rtype: bool
    """
    i = ord(ch)
    # A/a/C/c or G/g or T/t/U/u
    if (i | 0b1100010 == 0b1100011) or (i | 0b1100000 == 0b1100111) or (i | 0b1100001 == 0b1110101):
        return True
    else:
        return False


class IntToNucl(object):
    def __init__(self):
        self._map = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}

    def get(self, i):
        """
            Map: {0, 1, 2, 3} -> {A, C, T, G}
            @type i: int
            @rtype: str
        """
        if 0 <= i < 4:
            return self._map[i]
        else:
            raise TypeError('Value "%s" not allowed' % i)


def getBaseProb(qs):
    """
        Given a quality score, get the probability that the base is correct
        and the probability that different base X is correct.
        Map: qs -> P(base), P(other base X)

        @param qs: quality score
        @type qs: str
        @rtype: (float, float)
        @return: (pBaseCorrect, pOtherBaseX)
    """
    pBaseError = np.power(10, ((-1.) * (ord(qs) - 33)) / 10.)
    return 1. - pBaseError, pBaseError / 3.


class QSArray(object):
    def __init__(self, dna=None, qsArrayFq=None, qsA1=None, qsA2=None, pos1Within2=None, lenEmptyQSA=None):
        """
            QS Array represents a DNA sequence, handles operations on the underlying quality scores.
            A DNA sequence is represented by a distribution of (A,T,G,C) and a coverage for each position.

            Possible argument combinations
            - To represent an empty QS Array of length (lenEmptyQSA)
            - To represent a read (dna, qsArray)
            - To represent a consensus QS Array of two arrays (qsA1, qsA2, pos1Within2)

            @param dna: DNA sequence
            @param qsArrayFq: QS array as a FASTQ string
            @param pos1Within2: position of qsA1 within qsA2 (can be negative)
            @param lenEmptyQSA: length of an empty QSArray

            @type dna: str
            @type qsArrayFq: str
            @type qsA1: QSArray
            @type qsA2: QSArray
            @type pos1Within2: int
            @type lenEmptyQSA: int
        """
        # map: nucl -> int
        self._intToNucl = IntToNucl()

        # find out which constructor to use
        if qsA1 is not None and qsA2 is not None and pos1Within2 is not None:
            assert dna is None and qsArrayFq is None and lenEmptyQSA is None
            self._init1(qsA1, qsA2, pos1Within2)
        elif dna is not None and qsArrayFq is not None:
            assert qsA1 is None and qsA2 is None and pos1Within2 is None and lenEmptyQSA is None
            self._init2(dna, qsArrayFq)
        elif lenEmptyQSA is not None:
            assert dna is None and qsArrayFq is None and qsA1 is None and qsA2 is None and pos1Within2 is None
            # initialize an empty instance
            self._len = lenEmptyQSA
            self._probMatrix = self._getProbMatrix(self._len)
        else:
            raise TypeError('Wrong argument combination')

    def _getProbMatrix(self, columns):
        """
            @type columns: int
            @rtype: ndarray
        """
        return np.zeros((PROB_MATRIX_ROWS, columns), dtype=np.float64)

    def _init1(self, qs1, qs2, pos1Within2):
        """
            Init an instance representing a consensus array of two arrays.

            @type qs1: QSArray
            @type qs2: QSArray
            @type pos1Within2: int
        """
        # read R2 starts before R1
        if pos1Within2 < 0:
            tmp = qs1
            qs1 = qs2
            qs2 = tmp
            pos1Within2 = -pos1Within2

        # prob. matrix length
        self._len = max(qs2._len, qs1._len + pos1Within2)
        self._probMatrix = self._getProbMatrix(self._len)

        # before overlap, read from R2
        for r2 in range(pos1Within2):
            for i in range(PROB_MATRIX_ROWS):
                self._probMatrix[i][r2] = qs2._probMatrix[i][r2]
        r2 = pos1Within2

        # overlap
        overlapLen = min(qs1._len, qs2._len - pos1Within2)
        for r1 in range(overlapLen):
            for i in range(PROB_MATRIX_ROWS):
                self._probMatrix[i][r2] = qs1._probMatrix[i][r1] + qs2._probMatrix[i][r2]
            r2 += 1

        # after overlap
        if overlapLen == qs1._len:
            # read from R2
            for r2 in range(r2, r2 + qs2._len - pos1Within2 - overlapLen):
                for i in range(PROB_MATRIX_ROWS):
                    self._probMatrix[i][r2] = qs2._probMatrix[i][r2]

        else:
            # read from R1
            for r1 in range(overlapLen, qs1._len):
                for i in range(PROB_MATRIX_ROWS):
                    self._probMatrix[i][r2] = qs1._probMatrix[i][r1]
                r2 += 1

    def _init2(self, dna, qsArray):
        """
            Init an instance given a FASTQ entry (DNA seq and QS as string).

            @type dna: str
            @type qsArray: str
        """
        # QS array length
        m = len(qsArray)
        assert len(dna) == m
        self._len = m

        # A, C, T, G
        a = nuclToInt('A')
        c = nuclToInt('C')
        t = nuclToInt('T')
        g = nuclToInt('G')

        # rows: {A, C, T, G, coverage} col: {0, ..., m-1}
        self._probMatrix = self._getProbMatrix(m)

        # fill in the probability matrix
        for i, char, qs in zip(range(m), list(dna), list(qsArray)):

            # get base probabilities
            pBaseCorrect, pOtherBaseX = getBaseProb(qs)

            try:
                # base char int representation
                charInt = nuclToInt(char)

                # probabilities of individual bases
                for j in range(PROB_MATRIX_ROWS - 1):
                    if j == charInt:
                        self._probMatrix[j][i] = pBaseCorrect
                    else:
                        self._probMatrix[j][i] = pOtherBaseX

            # the base is ambiguous (a non-DNA characters)
            except TypeError:
                x = self._probMatrix
                ch = str(char).upper()
                pc = pBaseCorrect
                pe = 1. - pBaseCorrect

                if ch == 'K':
                    x[g][i] = x[t][i] = pc / 2.
                    x[c][i] = x[a][i] = pe / 2.
                elif ch == 'M':
                    x[a][i] = x[c][i] = pc / 2.
                    x[t][i] = x[g][i] = pe / 2.
                elif ch == 'R':
                    x[a][i] = x[g][i] = pc / 2.
                    x[t][i] = x[c][i] = pe / 2.
                elif ch == 'Y':
                    x[c][i] = x[t][i] = pc / 2.
                    x[g][i] = x[a][i] = pe / 2.
                elif ch == 'S':
                    x[c][i] = x[g][i] = pc / 2.
                    x[a][i] = x[t][i] = pe / 2.
                elif ch == 'W':
                    x[a][i] = x[t][i] = pc / 2.
                    x[g][i] = x[c][i] = pe / 2.
                elif ch == 'B':
                    x[c][i] = x[g][i] = x[t][i] = pc / 3.
                    x[a][i] = pe
                elif ch == 'V':
                    x[a][i] = x[c][i] = x[g][i] = pc / 3.
                    x[t][i] = pe
                elif ch == 'H':
                    x[a][i] = x[c][i] = x[t][i] = pc / 3.
                    x[g][i] = pe
                elif ch == 'D':
                    x[a][i] = x[g][i] = x[t][i] = pc / 3.
                    x[c][i] = pe
                else:
                    for j in range(PROB_MATRIX_ROWS - 1):
                        self._probMatrix[j][i] = 1. / (PROB_MATRIX_ROWS - 1)

            # set the coverage of the base
            self._probMatrix[PROB_MATRIX_ROWS - 1][i] = 1.

    def getLen(self):
        return self._len

    def getConsDna(self):
        """
            Get consensus DNA sequence based on the most probable bases.
            @rtype: str
            @return: the representative consensus DNA
        """
        dnaA = np.zeros(self._len, dtype=np.int8)
        for i in range(self._len):
            maxI = 0
            for j in range(1, PROB_MATRIX_ROWS - 1):
                if self._probMatrix[j][i] > self._probMatrix[maxI][i]:
                    maxI = j
            dnaA[i] = maxI

        return ''.join(map(lambda x: self._intToNucl.get(x), dnaA))

    def getQSStr(self, dna):
        """
            Get the QS FASTQ string for the representative DNA.
            @param dna: representative DNA of this QS array.
            @type dna: str
            @rtype: str
            @return: consensus QS for the representative DNA
        """
        qsA = np.zeros(self._len, dtype=np.int8)
        for i, c in zip(range(self._len), list(dna)):

            p = self._probMatrix[nuclToInt(c)][i] / self._probMatrix[PROB_MATRIX_ROWS - 1][i]
            qsA[i] = int(round((-10.) * np.log10(1. - p)))

        return ''.join(map(lambda x: chr(x + 33), qsA))

    def revCompl(self):
        """
            Make reverse complement of the QS array.
        """
        assert PROB_MATRIX_ROWS == 5

        i = 0
        j = self._len - 1

        tmp = np.zeros(5, dtype=np.float64)
        interval1 = [0, 1, 2, 3, 4]
        interval2 = [2, 3, 0, 1, 4]

        while i <= j:

            for k in range(5):
                tmp[k] = self._probMatrix[k][i]

            if i < j:
                for k, r in zip(interval1, interval2):
                    self._probMatrix[k][i] = self._probMatrix[r][j]

            for k, r in zip(interval1, interval2):
                self._probMatrix[k][j] = tmp[r]

            i += 1
            j -= 1

    def getPosCovArray(self):
        """
            @return: coverage of each base in the probability matrix
            @rtype: ndarray
        """
        return self._probMatrix[PROB_MATRIX_ROWS - 1]

    def getAvgCov(self):
        """
            @return: average coverage of the consensus sequence represented by the probability matrix
            @rtype: float
        """
        return self._probMatrix[PROB_MATRIX_ROWS - 1].sum() / float(len(self._probMatrix[PROB_MATRIX_ROWS - 1]))

    def getOverlapScore(self, selfPos, qsArray2, pos2, overlapLen):
        """
            Get the overlap score of two QS arrays considering matching DNA sequences.

            @attention: for the normalization, divide by the overlap-length.
            @type selfPos: int
            @type qsArray2: QSArray
            @type pos2: int
            @type overlapLen: int
            @rtype: (float, float)
        """
        probMatrix2 = qsArray2._probMatrix
        assert 0 <= selfPos < self._len and selfPos + overlapLen <= self._len
        assert 0 <= pos2 < qsArray2._len and pos2 + overlapLen <= qsArray2._len
        sumA = np.zeros(3, dtype=np.float128)
        for i in range(overlapLen):
            sumA[0] = 0.
            for j in range(PROB_MATRIX_ROWS - 1):
                covProduct = (self._probMatrix[PROB_MATRIX_ROWS - 1][selfPos] * probMatrix2[PROB_MATRIX_ROWS - 1][pos2])
                sumA[0] += (self._probMatrix[j][selfPos] * probMatrix2[j][pos2]) / covProduct

            sumA[1] += np.log10(sumA[0])
            sumA[2] += sumA[0]
            selfPos += 1
            pos2 += 1

        return (np.power(10, float(sumA[1] / overlapLen)), sumA[2])

    def getOverlapScoreProt(self, selfPos, qsArray2, pos2, overlapLen):
        """
            Get the overlap score of two QS arrays considering matching of PROT sequences.

            @attention: for the normalization, divide by the overlap-length.
            @type selfPos: int
            @type qsArray2: QSArray
            @type pos2: int
            @type overlapLen: int
            @rtype: (float, float)
        """
        assert overlapLen % 3 == 0
        assert PROB_MATRIX_ROWS >= 5
        assert 0 <= selfPos < self._len and selfPos + overlapLen <= self._len
        assert 0 <= pos2 < qsArray2._len and pos2 + overlapLen <= qsArray2._len

        sum = np.zeros(3, dtype=np.float128)
        a = nuclToInt('A')
        c = nuclToInt('C')
        t = nuclToInt('T')
        g = nuclToInt('G')

        x = self._probMatrix
        y = qsArray2._probMatrix
        s = PROB_MATRIX_ROWS - 1

        # sum up over all aminoacids
        for e in range(0, overlapLen, 3):

            i = e + selfPos
            j = e + pos2

            # K
            sum[0] = x[a][i] * x[a][i+1] * (x[a][i+2] + x[g][i+2]) \
                   * y[a][j] * y[a][j+1] * (y[a][j+2] + y[g][j+2])
            # N
            sum[0] += x[a][i] * x[a][i+1] * (x[c][i+2] + x[t][i+2]) \
                    * y[a][j] * y[a][j+1] * (y[c][j+2] + y[t][j+2])
            # T
            sum[0] += x[a][i] * x[c][i+1] * x[s][i+2] \
                    * y[a][j] * y[c][j+1] * y[s][j+2]
            # R
            sum[0] += (x[a][i] * x[g][i+1] * (x[a][i+2] + x[g][i+2]) + x[c][i] * x[g][i+1] * x[s][i+2]) \
                    * (y[a][j] * y[g][j+1] * (y[a][j+2] + y[g][j+2]) + y[c][j] * y[g][j+1] * y[s][j+2])
            # S
            sum[0] += (x[a][i] * x[g][i+1] * (x[c][i+2] + x[t][i+2]) + x[t][i] * x[c][i+1] * x[s][i+2]) \
                    * (y[a][j] * y[g][j+1] * (y[c][j+2] + y[t][j+2]) + y[t][j] * y[c][j+1] * y[s][j+2])
            # I
            sum[0] += x[a][i] * x[t][i+1] * (x[a][i+2] + x[c][i+2] + x[t][i+2]) \
                    * y[a][j] * y[t][j+1] * (y[a][j+2] + y[c][j+2] + y[t][j+2])
            # M
            sum[0] += x[a][i] * x[t][i+1] * x[g][i+2] \
                    * y[a][j] * y[t][j+1] * y[g][j+2]
            # H
            sum[0] += x[c][i] * x[a][i+1] * (x[c][i+2] + x[t][i+2]) \
                    * y[c][j] * y[a][j+1] * (y[c][j+2] + y[t][j+2])
            # Q
            sum[0] += x[c][i] * x[a][i+1] * (x[a][i+2] + x[g][i+2]) \
                    * y[c][j] * y[a][j+1] * (y[a][j+2] + y[g][j+2])
            # P
            sum[0] += x[c][i] * x[c][i+1] * x[s][i+2] \
                    * y[c][j] * y[c][j+1] * y[s][j+2]
            # L
            sum[0] += (x[c][i] * x[t][i+1] * x[s][i+2] + x[t][i] * x[t][i+1] * (x[a][i+2] + x[g][i+2])) \
                    * (y[c][j] * y[t][j+1] * y[s][j+2] + y[t][j] * y[t][j+1] * (y[a][j+2] + y[g][j+2]))
            # E
            sum[0] += x[g][i] * x[a][i+1] * (x[a][i+2] + x[g][i+2]) \
                    * y[g][j] * y[a][j+1] * (y[a][j+2] + y[g][j+2])
            # D
            sum[0] += x[g][i] * x[a][i+1] * (x[c][i+2] + x[t][i+2]) \
                    * y[g][j] * y[a][j+1] * (y[c][j+2] + y[t][j+2])
            # A
            sum[0] += x[g][i] * x[c][i+1] * x[s][i+2] \
                    * y[g][j] * y[c][j+1] * y[s][j+2]
            # G
            sum[0] += x[g][i] * x[g][i+1] * x[s][i+2] \
                    * y[g][j] * y[g][j+1] * y[s][j+2]
            # V
            sum[0] += x[g][i] * x[t][i+1] * x[s][i+2] \
                    * y[g][j] * y[t][j+1] * y[s][j+2]
            # *
            sum[0] += (x[t][i] * (x[a][i+1] * (x[a][i+2] + x[g][i+2]) + x[g][i+1] * x[a][i+2])) \
                    * (y[t][j] * (y[a][j+1] * (y[a][j+2] + y[g][j+2]) + y[g][j+1] * y[a][j+2]))
            # Y
            sum[0] += x[t][i] * x[a][i+1] * (x[t][i+2] + x[c][i+2]) \
                    * y[t][j] * y[a][j+1] * (y[t][j+2] + y[c][j+2])
            # C
            sum[0] += x[t][i] * x[g][i+1] * (x[t][i+2] + x[c][i+2]) \
                    * y[t][j] * y[g][j+1] * (y[t][j+2] + y[c][j+2])
            # W
            sum[0] += x[t][i] * x[g][i+1] * x[g][i+2] \
                    * y[t][j] * y[g][j+1] * y[g][j+2]
            # F
            sum[0] += x[t][i] * x[t][i+1] * (x[t][i+2] + x[c][i+2]) \
                    * y[t][j] * y[t][j+1] * (y[t][j+2] + y[c][j+2])

            # division ~ product of coverage
            d = x[s][i] * x[s][i+1] * x[s][i+2] * y[s][j] * y[s][j+1] * y[s][j+2]

            entry = sum[0] / d
            sum[1] += np.log10(entry)
            sum[2] += entry

        # other interesting statistics:
        # float(sum[1])
        # float(sum[1] / (overlapLen / 3))
        # (sum[2] / (overlapLen / 3))
        return (np.power(10, float(sum[1] / (overlapLen / 3))), sum[2] * 3)

    def update(self, startPos, qs):
        """
            Updates this QSArray from startPos in this array by another "qs" array.
            @type startPos:int
            @type qs: QSArray
        """
        assert 0 <= startPos < self._len and qs._len <= self._len - startPos
        for i in range(qs._len):
            for j in range(PROB_MATRIX_ROWS):
                self._probMatrix[j][startPos + i] += qs._probMatrix[j][i]

# TESTS ---------------------------------------------------

def getRandQS(length, rand=None):
    """
        @return: random quality score
    """
    if rand is None:
        rand = np.random.RandomState(length)
    l = []
    for i in range(length):
        l.append(chr(rand.randint(37) + 33 + 5))
    return ''.join(l)


def _test1():
    rand = np.random.RandomState(0)
    dom = 'GAAAGTTGACCAACTGATATTTGCCGGTCTGGCATCAAGTTATTCGGTATTGAGGGAAGATGAACGTGAACTGGGTGTCTGCGTCGTCGATATCGGTGGTGGTA' \
          'CAATGAATATCGCCGTTTATACCGGTGGGGCATTGCGCCACACTAAGGTAATTCCTTATGCTGGCAATGTAGTGACCAG'
    rare = 'GAAAGTTGACCAACTGATATcTGCCGGTCTGGCATCAAGTTATTCGGTAcTGAGGGAAGATGAACGTGAACTGGGTGTCTGCGTCGTCGATATCGGTGGTGGT' \
           'ACAATGAATATCGCCGTTTATACCtGTGGGGCATTGCGCCACACTAAGGTAATTCCTTATGCTGGCAATGTAGTGACCAG'
    length = len(dom)

    qsAll = QSArray(lenEmptyQSA=length)

    for i in range(50):
        qsAll.update(0, QSArray(dna=dom, qsArrayFq=getRandQS(length, rand)))

    for i in range(10):
        qsAll.update(0, QSArray(dna=rare, qsArrayFq=getRandQS(length, rand)))

    print 'prefix: dom vs align.'
    qs2 = QSArray(dna=dom[:60+57], qsArrayFq=getRandQS(60+57, rand))
    print qsAll.getOverlapScore(0, qs2, 0, 60+57)

    print 'prefix: rare vs align.'
    qs2 = QSArray(dna=rare[:60+57], qsArrayFq=getRandQS(60+57, rand))
    print qsAll.getOverlapScore(0, qs2, 0, 60+57)

    print 'suffix: dom vs align.'
    qs2 = QSArray(dna=dom[60:], qsArrayFq=getRandQS(57+66, rand))
    print qsAll.getOverlapScore(60, qs2, 0, 57+66)

    print 'suffix: rare vs align.'
    qs2 = QSArray(dna=rare[60:], qsArrayFq=getRandQS(57+66, rand))
    print qsAll.getOverlapScore(60, qs2, 0, 57+66)


def _test2():
    dna = 'ATGCTACGAACACAGCAACCA'
    qs = 'G(C#$1GC8GCG=GG=1CG=C'
    qsA = QSArray(dna=dna, qsArrayFq=qs)

    rc = str(Seq(dna, generic_dna).reverse_complement())
    print '(+)strand prob. matrix:\n%s' % qsA._probMatrix
    qsA.revCompl()
    print '(-)strand prob. matrix:\n%s' % qsA._probMatrix
    print '(+)strand\n%s\n(-)strand\n%s\n%s' % (dna, rc, qsA.getConsDna())
    assert rc == qsA.getConsDna()
