"""
    FASTA read/write module version: 1.2

    The MIT License (MIT)

    Copyright (c) 2014  Ivan Gregor

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
    documentation files (the "Software"), to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions
    of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
    THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    Contains basic functionality to work with FASTA files.
"""

import sys
import os
import re
import types
import gzip
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from algbioi.com.csv import OutFileBuffer
from algbioi.com.common import removeNonDna
from algbioi.com.common import noNewLine


def dnaToProt(inFastaDna, outFastaProt, translTable=11):
    """
        Translates DNA gene! sequences to PROT sequences.

        @param inFastaDna: input fasta file containing DNA sequences
        @param outFastaProt: output fasta file containing sequences translated to protein sequences
        @param translTable: default 11 for bacteria and archaea (http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    """
    out = OutFileBuffer(outFastaProt)
    for seqName, seqDna in fastaFileToDictWholeNames(inFastaDna).iteritems():
        seqProt = Seq(seqDna, generic_dna).translate(table=translTable, stop_symbol='', cds=True)
        out.writeText('>%s\n%s\n' % (seqName, seqProt))
    out.close()


def fastaToProt(inFasta, outFasta, translTable=11):
    """
        Translate sequences from the input fasta file to all six reading frames and store them into a FASTA file.
    """
    out = OutFileBuffer(outFasta)

    for name, seq in fastaFileToDictWholeNames(inFasta).iteritems():
        sLen = len(seq)
        revSeq = str(Seq(seq, generic_dna).reverse_complement())

        assert len(seq) == len(revSeq)

        out.writeText('>%s_1\n%s\n' % (name, Seq(seq[:(sLen / 3) * 3], generic_dna).translate(table=translTable)))
        out.writeText('>%s_2\n%s\n' % (name, Seq(seq[1:((sLen - 1) / 3) * 3 + 1], generic_dna).translate(table=translTable)))
        out.writeText('>%s_3\n%s\n' % (name, Seq(seq[2:((sLen - 2) / 3) * 3 + 2], generic_dna).translate(table=translTable)))

        out.writeText('>%s_4\n%s\n' % (name, Seq(revSeq[:(sLen / 3) * 3], generic_dna).translate(table=translTable)))
        out.writeText('>%s_5\n%s\n' % (name, Seq(revSeq[1:((sLen - 1) / 3) * 3 + 1], generic_dna).translate(table=translTable)))
        out.writeText('>%s_6\n%s\n' % (name, Seq(revSeq[2:((sLen - 2) / 3) * 3 + 2], generic_dna).translate(table=translTable)))

    out.close()


def getSeqFromDir(seqId, dirList, baseFileName):
    """
        Search for a sequence (seqId) in all files defined by a baseFileName and a list of directories,
        return the first sequence found.

        @param seqId: sequence id of the sequence contained in one of the files defined by (dirList, baseFileName)
        @param dirList: list of directories in which it searches for the fasta file containing the seqId
        @param baseFileName: defines the fasta files for the search (directory/baseFileName.fna[.gz])
        @type seqId: str
        @type dirList: list[str]
        @type baseFileName: str
        @return: a sequence or None
        @rtype: str
    """
    for d in dirList:
        assert os.path.isdir(d)
        filePath = os.path.join(d, baseFileName + '.fna')
        if not os.path.isfile(filePath):
            filePath += '.gz'
            if not os.path.isfile(filePath):
                filePath = None
        if filePath is not None:
            seqIdToSeq = fastaFileToDictWholeNames(filePath)
            if seqId in seqIdToSeq:
                return seqIdToSeq[seqId]
    return None


def getSequenceBuffer(dirList):
    """
        Reads all fasta files in the list of directories and stores the mapping in a dictionary (thread safe).
        Here a fasta file name (without suffix ".fna" or ".fna.gz") is interpreted as a strain id.

        @param dirList: list of directories containing fasta files
        @type dirList: list[str]
        @rtype: dict[(str,str),str]
        @return: map: (strainId, seqId) -> sequence
    """
    man = multiprocessing.Manager()
    mapping = man.dict()
    for d in dirList:
        assert os.path.isdir(d)
        for f in os.listdir(d):
            path = os.path.join(d, f)
            if os.path.isfile(path) and (path.endswith('.fna') or path.endswith('.fna.gz')):
                if path.endswith('.fna'):
                    strainId = f[:-4]
                else:
                    strainId = f[:-7]
                for seqName, seq in fastaFileToDictWholeNames(path).iteritems():
                    k = (strainId, seqName)
                    assert k not in mapping
                    mapping[k] = seq
    return mapping


def sortSeqDesc(inFasta, outFasta):
    """
        Sort sequences in a descending order.

        @param inFasta: input fasta file
        @param outFasta: output sorted fasta file
    """
    tupleList = []
    for seqName, seq in getSequencesToList(inFasta):
        tupleList.append((seqName, seq, len(seq)))
    # sort
    tupleList.sort(key=lambda x: x[2], reverse=True)

    out = OutFileBuffer(outFasta)
    for seqName, seq, bp in tupleList:
        out.writeText('>%s\n%s\n' % (seqName, seq))
    out.close()


def cmpSeqFiles(filePath1, filePath2, verbose=False, format='fastq'):
    """
        Compares two sequence files.

        @attention: uses SeqIO.parse, thus can be slow for very large files

        @return: True if both files contain the same entries, else False.
    """
    d1 = {}
    d2 = {}
    f1 = open(filePath1)
    f2 = open(filePath2)
    for record in SeqIO.parse(f1, format):
        d1[record.id] = record
    for record in SeqIO.parse(f2, format):
        d2[record.id] = record
    f1.close()
    f2.close()
    if len(d1) != len(d2):
        if verbose:
            print('Different lengths %s %s' % (len(d1), len(d2)))
        return False
    for k, v1 in d1.iteritems():
        v2 = d2[k]
        if str(v1) != str(v2) \
                or str(v1.letter_annotations['phred_quality']) != str(v2.letter_annotations['phred_quality']):
            if verbose:
                print('Different sequences! %s %s' % (v1, v2))
            return False
    if verbose:
        print('Files contain the same sequences: %s %s' % (filePath1, filePath2))
    return True


def splitPairedReads(inPairedFasta, outEvenFasta, outOddFasta):
    _forEachRecord(inPairedFasta, SplitFasta(outEvenFasta, outOddFasta)).close()


class SplitFasta():
    def __init__(self, evenFasta, oddFasta):
        self._evenFasta = OutFileBuffer(evenFasta)
        self._oddFasta = OutFileBuffer(oddFasta)
        self._counter = 0

    def parse(self, record):
        entry = '>' + str(record.id) + '\n' + str(record.seq) + '\n'
        if self._counter % 2 == 0:
            self._evenFasta.writeText(entry)
        else:
            self._oddFasta.writeText(entry)
        self._counter += 1

    def close(self):
        self._oddFasta.close()
        self._evenFasta.close()


def filterOutSequences(inFileName, outFileName, allowedNamesSet, formatName="fasta", seqNameModifyFunction=None):
    """
        From the input fasta file filter out sequences their names are not contained in the allowedNamesSet.

        @param allowedNamesSet: the set of entries that are allowed as a sequence names
        @param seqNameModifyFunction: a sequence`s name is modified by this function and then compared to the
        allowedNamesSet
    """
    outFileBuffer = OutFileBuffer(outFileName)
    recordCondition = RecordConditionFilterOutSequences(allowedNamesSet, seqNameModifyFunction)
    parser = RecordFilter(outFileBuffer, formatName, recordCondition)
    _forEachRecord(inFileName, parser)


def filterOutNonDna(inFileName, outFileName):
    outFileBuffer = OutFileBuffer(outFileName)
    parser = RemoveNonDnaParser(outFileBuffer)
    _forEachRecord(inFileName, parser)


def getSequenceToBpDict(fastaFilePath):
    """
        Reads a fasta file and returns mapping: sequenceName -> sequenceLength.
    """
    return _forEachRecord(fastaFilePath, SeqToBpParser()).getSeqToBpDict()


class SeqToBpParser():
    def __init__(self):
        self._seqToBp = dict([])

    def parse(self, record):
        self._seqToBp[record.id] = len(str(record.seq))

    def getSeqToBpDict(self):
        return self._seqToBp

    def getFormatName(self):
        return "fasta"


def getSequencesToList(fastaFilePath):
    """
        Reads a fasta file and returns a list of: (sequenceName, sequence).
    """
    return _forEachRecord(fastaFilePath, SeqToListParser()).getSeqToList()


class SeqToListParser():
    def __init__(self):
        self._seqToList = []

    def parse(self, record):
        self._seqToList.append((str(record.id), noNewLine(str(record.seq))))

    def getSeqToList(self):
        return self._seqToList

    def getFormatName(self):
        return "fasta"


class RemoveNonDnaParser():
    def __init__(self, outFileBuffer):
        self._outFileBuffer = outFileBuffer

    def finalize(self):
        self._outFileBuffer.close()

    def getFormatName(self):
        return "fasta"

    def parse(self, record):
        self._outFileBuffer.writeText(str('>' + str(record.id) + '\n'))
        self._outFileBuffer.writeText(str(removeNonDna(str(record.seq)) + '\n'))


class RecordConditionFilterOutSequences():
    def __init__(self, allowedNamesSet, seqNameModifyFunction=None):
        self.allowedNamesSet = allowedNamesSet
        self.seqNameModifyFunction = seqNameModifyFunction

    def takeRecord(self, record):
        """
            If the record.id (modified by the function) is in the allowedNamesSet then the entry will be accepted.
        """
        idr = record.id
        if self.seqNameModifyFunction is not None:
            idr = self.seqNameModifyFunction(idr)
        if idr in self.allowedNamesSet:
            return True
        else:
            return False


class RecordFilter():
    """
        Appends a record that is currently parsed to the outFileBuffer if it satisfies the condition.
    """
    def __init__(self, outFileBuffer, formatName, recordCondition):
        self.outFileBuffer = outFileBuffer
        self.formatName = formatName
        self.recordCondition = recordCondition

    def parse(self, record):
        if self.recordCondition.takeRecord(record):
            self.outFileBuffer.writeText(record.format(self.formatName))

    def getFormatName(self):
        return self.formatName

    def finalize(self):
        self.outFileBuffer.close()


def fastaFileToDict(fastaFilePath, formatName='fasta'):
    """
        Reads a fasta file and returns mapping: seqName -> sequence.
    """
    return _forEachRecord(fastaFilePath, _RecordStorage(formatName), formatName=formatName).getSeqNameToSeq()


def cpSeqNoShortSeq(inFile, outFile, minLen):
    """
        Copy sequences longer or equal to a minimum length from the input to the output file.

        @param inFile: input fasta file
        @param outFile: output fasta file containing only sequences longer or equal to the minimum length
        @param minLen: minimum length of a sequence that will be copied
    """
    out = OutFileBuffer(outFile)
    first = True
    for name, seq in fastaFileToDictWholeNames(inFile).iteritems():
        if len(seq) >= minLen:
            if first:
                out.writeText('>%s\n%s' % (name, seq))
                first = False
            else:
                out.writeText('\n>%s\n%s' % (name, seq))
    out.close()


def fastaFileToDictWholeNames(filePath):
    """
        Reads a fasta file and returns mapping: seqName -> sequence the whole sequence name is used
        as seqName!!! (even if it contains space)
    """
    seqIdToSeq = {}
    f = None
    try:
        if filePath.endswith('.gz'):
            f = gzip.open(os.path.normpath(filePath), mode='r')
        else:
            f = open(os.path.normpath(filePath), 'r')
    except Exception:
        print "Cannot open file:", filePath
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqIdToSeq[name] = seq
                    seq = ''
                name = line.replace('>', '')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqIdToSeq[name] = seq
    finally:
        if f is not None:
            f.close()
    return seqIdToSeq


class _RecordStorage():
    def __init__(self, formatName='fasta'):
        self._seqNameToSeq = {}
        self._formatName = formatName

    def parse(self, record):
        self._seqNameToSeq[str(record.id)] = str(record.seq)

    def getFormatName(self):
        return self._formatName

    def getSeqNameToSeq(self):
        return self._seqNameToSeq


def _forEachRecord(filePath, parser, formatName="fasta"):
    """
        Call the parser for each record in the file.
    """
    try:
        if isinstance(parser.getFormatName, types.MethodType):
            formatName = parser.getFormatName()
    except Exception:
        pass
    try:
        f = open(os.path.normpath(filePath), 'r')
    except Exception:
        sys.stderr.write('Cannot open a %s file for reading: %s\n' % (formatName, filePath))
        raise
    else:
        try:
            readBuffer = SeqIO.parse(f, parser.getFormatName())
            for record in readBuffer:
                parser.parse(record)
        except Exception:
            sys.stderr.write('Cannot read from a ' + formatName + ' file: ' + filePath + '\n')
            raise
        finally:
            f.close()
    try:
        if isinstance(parser.finalize, types.MethodType):
            parser.finalize()
    except Exception:
        pass

    return parser
