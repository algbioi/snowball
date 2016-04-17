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

    Represents a (super-)read record (i.e. a record representing a consensus sequence).
    Module version 1.2
"""

import copy
from algbioi.haplo import qs as qs_man


class ReadRec(object):
    def __init__(self, recordId, qsArray, readingFrameTag, annotStart, annotLen, hmmCoordStart, hmmCoordLen=None,
                 protSeq=None, protStart=None, protLen=None):
        """
            A record representing one read (also a super-read or a joined pair-end read).
            The constructor is used to create an initial record, the function "merge records" of this module
            is then used to create super-read-records by merging of read-records.

            @param recordId: record identifier
            @param qsArray: QS array representing a DNA sequence of the (super-)read (oriented as +strand)
            @param readingFrameTag: reading frame indicator (1..6), (for super-reads: 1..3)

            @param annotStart: start of the HMM annotation within the read (0 based)
            @param annotLen: length of the HMM annotation
            @param hmmCoordStart: start of the HMM annotation within the Gene family (0 based)
            @param hmmCoordLen: length of the annotation withing the Gene family (or None)

            @param protSeq: protein sequence (or None)
            @param protStart: start of the HMM annotation on the PROT sequence (0 based, or None)
            @param protLen: length of the HMM annotation on the PROT sequence (or None)

            @type recordId: str
            @type qsArray: qs_man.QSArray
            @type readingFrameTag: int
            @type annotStart: int
            @type annotLen: int
            @type hmmCoordStart: int
            @type hmmCoordLen: int | None
            @type protSeq: str | None
            @type protStart: int | None
            @type protLen: int | None
        """
        assert annotLen == 3 * protLen

        # if the gene family annotation is on the reverse strand, orient it as if it was on the positive strand
        if 4 <= readingFrameTag <= 6:
            # dna reverse complement
            qsArray.revCompl()
            # reverse the start annotation position
            annotStart = qsArray.getLen() - annotStart - annotLen

        self.recordId = recordId
        self.qsArray = qsArray
        self.readingFrameTag = readingFrameTag
        self.annotStart = annotStart
        self.annotLen = annotLen
        self.hmmCoordStart = hmmCoordStart
        self.hmmCoordLen = hmmCoordLen
        self.protSeq = protSeq
        self.protStart = protStart
        self.protLen = protLen

        # DNA sequence of the record
        self.dnaSeq = qsArray.getConsDna()

        # a list of ReadRec merging of which this read-record originated
        # None - if this read-record did not originate by merging of other read-records
        # (this is used only for the evaluation, can be removed to speedup the algorithm)
        self.readRecList = None

        # the position of this (super-)read within a bigger contig (super-read)
        self.posWithinContig = 0

        # label for the evaluation purposes (this is used only for the evaluation as well, i.e. can be removed)
        self.labelEval = None

        # hmm coordinates for evaluation purposes
        self.evalHmmCoord = None

    def getPosCovArray(self):
        """
            @rtype: ndarray
            @return: position coverage array, entries as floats! (i.e. for each position,
            how many pair-end reads cover it, e.g. 111222111 for a joined pair-end read)
        """
        return self.qsArray.getPosCovArray()

    def getLen(self):
        """
            @return: length of the (super-)read
            @rtype: int
        """
        return self.qsArray.getLen()

    def isSuperRead(self):
        """
            @rtype: bool
            @return: is this a super-read
        """
        return self.readRecList is not None

    def getConstituentReadCount(self):
        """
            @return: from how many read records this read record consist of
            @rtype: int
        """
        if self.readRecList is None:
            return 1
        else:
            return len(self.readRecList)

    def getAvgCov(self):
        """
            @return: (super-)read average coverage
        """
        return self.qsArray.getAvgCov()

    def getLabels(self):
        """
            For the evaluation.
            @rtype: list | None
        """
        try:
            return self.labelEval
        except AttributeError:
            return None


def mergeRecords(rec1, rec2, overlapIdx):
    """
        Merges record-1 into record-2, returns a merged object.

        @param rec1: (super-)read-record
        @param rec2: read-record (that hasn't been joined with any other so far)
        @param overlapIdx: position within rec2, start pos of rec1 mapped onto rec2 (it can be a negative value)
        @param score: score used for this overlap

        @type rec1: ReadRec
        @type rec2: ReadRec
        @type overlapIdx: int
        @type score: int
        @rtype: ReadRec

        @return: a merged read record
    """
    # make a shallow copy of the object
    if rec2.readRecList is None:
        mRec = copy.copy(rec1)
    else:
        mRec = copy.copy(rec2)

    if overlapIdx < 0:
        tmp = rec1
        rec1 = rec2
        rec2 = tmp
        overlapIdx = -overlapIdx

    # concatenate record ids (this can be simplified)
    mRec.recordId = rec1.recordId + '|' + rec2.recordId

    # part before overlap rec1
    boLen1 = 0

    # part before overlap rec2
    boLen2 = overlapIdx

    # get the consensus QS Array
    mRec.qsArray = qs_man.QSArray(qsA1=rec1.qsArray, qsA2=rec2.qsArray, pos1Within2=overlapIdx)

    # consensus DNA
    mRec.dnaSeq = mRec.qsArray.getConsDna()

    mRec.posWithinContig = 0

    # set the reading frame tag
    mRec.readingFrameTag = rec2.readingFrameTag

    if mRec.readingFrameTag > 3:
        mRec.readingFrameTag -= 3

    # Hmm annotation start positions
    anPos1 = rec1.annotStart + boLen2
    anPos2 = rec2.annotStart + boLen1

    # Hmm annotation end + 1 positions
    anPosEndA1 = anPos1 + rec1.annotLen
    anPosEndA2 = anPos2 + rec2.annotLen

    # the annotations must have an overlap
    assert max(anPos1, anPos2) < min(anPosEndA1, anPosEndA2)

    # get the union of the annotations
    mRec.annotStart = min(anPos1, anPos2)
    mRec.annotLen = max(anPosEndA1, anPosEndA2) - mRec.annotStart

    # Hmm coord within gene family, take coordinates of the one at the interval boundary or the one with highest support
    if anPos1 < anPos2:
        mRec.hmmCoordStart = rec1.hmmCoordStart
    elif anPos1 > anPos2:
        mRec.hmmCoordStart = rec2.hmmCoordStart
    elif rec1.annotLen >= rec2.annotLen:
        mRec.hmmCoordStart = rec1.hmmCoordStart
    else:
        mRec.hmmCoordStart = rec2.hmmCoordStart

    # a list of all read-records this merged read-record consist of
    if rec2.readRecList is None:

        if rec1.readRecList is None:

            # two reads were merged
            mRec.readRecList = [rec1, rec2]

            # set the positions of the read-records within the merged record
            rec1.posWithinContig = boLen2
        else:
            # super-read and a read were merged
            mRec.readRecList = rec1.readRecList + [rec2]

            # set the positions of the read-records within the merged record
            for r in rec1.readRecList:
                r.posWithinContig += boLen2

        rec2.posWithinContig = 0
    else:

        if rec1.readRecList is None:
            # read and super-read were merged
            mRec.readRecList = rec2.readRecList + [rec1]
            rec1.posWithinContig = boLen2
        else:
            # two super reads were merged
            mRec.readRecList = rec1.readRecList + rec2.readRecList
            for r in rec1.readRecList:
                r.posWithinContig += boLen2

    # these entries are currently not used for super-reads (may be dropped or used in the future)
    mRec.protSeq = None
    mRec.protStart = None
    mRec.protLen = None
    mRec.hmmCoordLen = None

    return mRec
