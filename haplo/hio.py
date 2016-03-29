#!/usr/bin/env python

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

    Parsing of the input files.
"""

import os
# import sys
import gzip
import cPickle
# import numpy as np

# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna

from algbioi.com import fq
from algbioi.com import fasta as fas
from algbioi.com import qs as qs_man
from algbioi.haplo import read_rec
from algbioi.hsim import pfam


def readDomblout(inDomblout):
    """
        Read in a dom file, allow only one hit per read.

        @return: map: readName -> (list of line tokens, frame-tag)
        @rtype: dict[str,(list[str], int)]
    """
    nameToHit = {}
    for line in gzip.open(inDomblout):
        line = line.strip()
        if line.startswith('#'):
            continue
        assert line.startswith('@')
        tokens = line.split()
        name = tokens[0][1:-2]
        tag = int(tokens[0][-1])

        assert name not in nameToHit
        nameToHit[name] = (tokens, tag)

    return nameToHit


def parse(inFq, inDomtblout, inProtFna=None):
    """
        Read in joined pair-end reads and its HMM annotation from the FASTQ, DOMTBLOUT, (and optionally PROT FASTA file)

        @param inFq: FASTQ file containing joined pair-end reads
        @param inDomtblout: HMM annotation file
        @param inProtFna: corresponding prot sequence (can be None)

        @type inFq: str
        @type inDomtblout: str
        @type inProtFna: str | None
        @rtype: list[read_rec.ReadRec]

        @return: a list of read-records
    """
    recList = []

    # read in prot sequences
    if inProtFna is None:
        nameToProtSeq = {}
    else:
        nameToProtSeq = fas.fastaFileToDictWholeNames(inProtFna)

    # read in dom file
    nameToDom = readDomblout(inDomtblout)

    assert inProtFna is None or len(nameToProtSeq) == len(nameToDom)

    # read in pair-end reads, create ReadRec
    for readName, dna, comment, qs in fq.ReadFqGen(inFq):

        readName = readName[1:]  # strip starting @

        protSeq = nameToProtSeq.get(readName, None)

        hit, frameTag = nameToDom[readName]

        annotStart, annotLen, strain, score, acc = pfam.dnaHitInfo(hit, dna, protSeq)

        if strain == 1:
            assert 1 <= frameTag <= 3
        else:
            assert 4 <= frameTag <= 6

        # TODO: consider also strain, score, acc ?

        # alignment env coord
        protStart = int(hit[19]) - 1
        protLen = int(hit[20]) - protStart

        assert annotLen == 3 * protLen

        # alignment coord
        protStartAli = int(hit[17]) - 1
        protLenAli = int(hit[18]) - protStartAli

        # hmm coordinates
        hmmCoordStart = int(hit[15]) - 1
        hmmCoordLen = int(hit[16]) - hmmCoordStart

        # the env coordinates start (and end) often before the alignment coordinates and end after, get the offsets
        offsetEnv = protStartAli - protStart
        offsetEnvE = protLen - offsetEnv - protLenAli

        assert offsetEnv >= 0
        assert offsetEnvE >= 0

        hmmCoordStart -= offsetEnv
        hmmCoordLen += offsetEnv + offsetEnvE

        # create the coverage array (NO NEED)
        # posCovArray = np.zeros(len(dna), dtype=np.uint8)
        # start = len(dna) - pairEndReadLen  # incl.
        # end = pairEndReadLen - 1  # incl.
        # for i in range(len(posCovArray)):
        #     if start <= i <= end:
        #         posCovArray[i] = 2
        #     else:
        #         posCovArray[i] = 1

        tokens = comment.split('\t')

        if len(tokens) == 5:
            # get the ends of the pair-end read and corresponding quality-scores
            p, dna1, qs1, dna2, qs2 = tokens
            # get the QSArray representation
            qsA1 = qs_man.QSArray(dna=dna1, qsArrayFq=qs1)

            # get the reverse complement of the second end of pair-end read
            qsA2 = qs_man.QSArray(dna=dna2, qsArrayFq=qs2)
            qsA2.revCompl()

            # dna2 = str(Seq(dna2, generic_dna).reverse_complement())
            # qs2 = qs2[::-1]
            # qsA2 = qs_man.QSArray(dna=dna2, qsArrayFq=qs2)

            # get the consensus QS Array representing of the joined read
            qsA = qs_man.QSArray(qsA1=qsA1, qsA2=qsA2, pos1Within2=len(dna2) - len(dna))
        else:
            # there is just a simple dna sequence (i.e. not joined reads)
            qsA = qs_man.QSArray(dna=dna, qsArrayFq=qs)

        # cDna = qsA.getConsDna()
        # the consensus sequences can differ in the case both read-ends have different char with the same quality-score
        # if cDna != dna:
        #     print dna
        #     print cDna
        #     print qs1
        #     print (' ' * (len(dna) - len(dna2))) + qs2[::-1]
        #     print qsA.getQSStr(cDna)
        #     print('')
        #     sys.exit(0)
        recList.append(read_rec.ReadRec(readName, qsA, frameTag, annotStart, annotLen, hmmCoordStart, hmmCoordLen,
                                        protSeq, protStart, protLen))

    return recList


def storeReadRec(recList, outFilePath, compressLevel=1):
    """
        Store a list of read-records to a file.

        @param recList: list of read-records to be stored
        @param outFilePath: output dmp file
        @type recList: list[read_rec.ReadRec]
        @type outFilePath: str
    """
    storeObj(recList, outFilePath, compressLevel)
    # open file for writing
    # out = gzip.open(outFilePath, 'wb', compressLevel)
    # write the list to the file
    # cPickle.dump(recList, out, cPickle.HIGHEST_PROTOCOL)
    # out.close()


def storeObj(obj, outFilePath, compressLevel=1):
    """
        Store an object to a file.
        @type obj: object
        @type outFilePath: str
    """
    # open file for writing
    out = gzip.open(outFilePath, 'wb', compressLevel)
    # write the object to the file
    cPickle.dump(obj, out, cPickle.HIGHEST_PROTOCOL)

    out.close()


def loadReadRec(srcFilePath):
    """
        Load a list of records from a file.

        @param srcFilePath: a file containing stored read-records
        @type srcFilePath: str
        @return: list of read-records or None
        @rtype: list[read_rec.ReadRec]
    """
    return loadObj(srcFilePath)
    # out = gzip.open(srcFilePath, 'rb')
    # try:
    #     return cPickle.load(out)
    # except EOFError:
    #     return None


def loadObj(srcFilePath):
    """
        Load an object from a file.

        @type srcFilePath: str
        @return: an object or None
        @rtype: object | None
    """
    out = gzip.open(srcFilePath, 'rb')
    try:
        return cPickle.load(out)
    except EOFError:
        return None

