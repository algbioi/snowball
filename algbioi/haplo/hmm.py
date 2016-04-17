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

    HMM results mapping. Module version 1.2
"""

import gzip

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


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


def dnaHitInfo(hit, dnaSeq, protSeq=None, translTable=11):
    """
        Given the hit entry from the domtblout file for a PROT sequence, return the hit coordinates on the DNA sequence.

        @type hit: list[str]
        @type dnaSeq: str

        @param protSeq: if equal to a sequence, verification is performed (no verification if None)

        @return: tuple (startOnRead, overlapLen, strain, score, acc)
    """
    # read hit info
    score = float(hit[13])
    fromEnv = int(hit[19]) - 1
    lenEnv = int(hit[20]) - fromEnv
    acc = float(hit[21])
    frame = int(hit[0][-1:])
    translSeq = None

    # length of the overlap on the dna sequence
    overlapLen = lenEnv * 3

    # frame shift on the positive strand
    if 1 <= frame <= 3:
        # start position on the DNA read
        startOnRead = fromEnv * 3 + frame - 1
        strain = 1

        # get the prot sequence based on the position on the DNA read
        if protSeq is not None:
            translSeq = str(Seq(dnaSeq[startOnRead:startOnRead + overlapLen], generic_dna).translate(translTable))

    else:
        assert 4 <= frame <= 6
        offset = frame - 4

        # start position on the DNA read
        startOnRead = len(dnaSeq) - 1 - offset - (fromEnv * 3 + overlapLen - 1)
        strain = -1

        # get the prot sequence based on the position on the DNA read
        if protSeq is not None:
            subStrRev = str(Seq(dnaSeq[startOnRead:(startOnRead + overlapLen)], generic_dna).reverse_complement())

            startOnRead2 = len(dnaSeq) - startOnRead - overlapLen
            subStrRev2 = str(Seq(dnaSeq, generic_dna).reverse_complement())[startOnRead2:(startOnRead2 + overlapLen)]
            assert subStrRev == subStrRev2

            translSeq = str(Seq(subStrRev, generic_dna).translate(translTable))

    # verify that the translated part of the DNA sequence really correspond to the substring of the prot sequence
    if protSeq is not None and translSeq != protSeq[fromEnv: fromEnv + lenEnv]:
        raise Exception('Translation not equal:\n%s\n%s\n%s' % (translSeq, protSeq[fromEnv: fromEnv + lenEnv], strain))

    return startOnRead, overlapLen, strain, score, acc
