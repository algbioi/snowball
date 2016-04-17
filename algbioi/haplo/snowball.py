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

    Implements the main "Snowball" algorithm of the package. Module version 1.2
"""

from algbioi.com import common as com
from algbioi.haplo import read_rec


class JoinInfo(object):
    def __init__(self, score, overlap, seedRec, neighborRec):
        """
            Represents a helper record representing an overlap of two (super-)reads.
        """
        self.score = score
        self.overlap = overlap
        self.seedRec = seedRec
        self.neighborRec = neighborRec


def getBestOverlap(inspectList, tripletMap, considerProtSeqCmp, pOverlapMin, scoreOverlapLenMin,
                   scoreOverlapAnnotLenMin, stopOverlapMaxMismatch, considerOnlyPOverlap=False):
    """
        Compute overlaps for all the record pairs and return the best one or None.

        @param inspectList: a list of tuples of records (seed, neighbor)
        @param tripletMap: map: triplet -> aminoacid
        @param considerProtSeqCmp: consider overlapping annotated regions of the sequences as protein sequences
        @param pOverlapMin: minimum overlap probability
        @param scoreOverlapLenMin: minimum length of the overlap (expected true length)
        @param scoreOverlapAnnotLenMin: minimum length of the annotated overlap (expected true length)
        @param stopOverlapMaxMismatch: stop inspecting an overlap if the fraction of mismatches exceeds this threshold

        @type inspectList: list[(read_rec.ReadRec, read_rec.ReadRec)]
        @type tripletMap: dict[str,str]
        @type considerProtSeqCmp: bool
        @type pOverlapMin: float
        @type scoreOverlapLenMin: float
        @type scoreOverlapAnnotLenMin: float
        @type stopOverlapMaxMismatch: float
        @rtype: JoinInfo | None

        @return: best matching records overlap or None
    """
    # list of read-record pairs that can be joined
    candidateJoinList = []

    # go over all read-record pairs and inspect its possible overlap scores
    for seedRec, neighborRec in inspectList:

        minRecLen = min(seedRec.getLen(), neighborRec.getLen())
        minScore = int(float(minRecLen) * scoreOverlapLenMin)
        minAnnotScore = int(float(minRecLen) * scoreOverlapAnnotLenMin)

        # find possible overlap positions
        overlapList = getPossibleOverlaps(seedRec, neighborRec)

        # there are possible overlap positions to evaluate
        if overlapList is not None:

            # inspect the overlap list, store possible overlaps
            possibleJoinList = []
            for overlap in overlapList:

                # ge the score of this overlap
                overlapProb, overlapScore, annotScore = inspectOverlap(seedRec, neighborRec, overlap, tripletMap,
                                                                        considerProtSeqCmp, minScore, minAnnotScore,
                                                                        stopOverlapMaxMismatch, continuousOverlap=True)
                # the scores are sufficient (good enough)
                if overlapProb is not None and overlapProb >= pOverlapMin \
                        and overlapScore >= minScore and annotScore >= minAnnotScore:

                    if considerOnlyPOverlap:
                        possibleJoinList.append((overlapProb, overlap))
                    else:
                        possibleJoinList.append((overlapProb * overlapScore, overlap))

            # there was at least one sufficient overlap of these two read-records
            if len(possibleJoinList) > 0:
                # get the overlap with the best score
                possibleJoinList.sort(key=lambda x: x[0], reverse=True)
                scoreB, overlapB = possibleJoinList[0]
                candidateJoinList.append(JoinInfo(scoreB, overlapB, seedRec, neighborRec))

    # return the read-record pair with the highest score
    if len(candidateJoinList) > 0:
        candidateJoinList.sort(key=lambda x: x.score, reverse=True)
        return candidateJoinList[0]
    else:
        # no sufficient overlap was found
        return None


def inspectOverlap(rec1, rec2, overlapIdx, tripletMap, considerProtSeqCmp, overlapLenMin, overlapAnnotLenMin,
                    stopOverlapMaxMismatch, continuousOverlap=False):
    """
        Compute the overlap probability and expected length of the matching positions in the overlap.

        @param rec1: first read-record
        @param rec2: second read-record
        @param overlapIdx: an overlap to be inspected, start pos of rec1 mapped onto rec2 (it can be a negative value)
        @param tripletMap: map: triplet -> aminoacid
        @param considerProtSeqCmp: consider overlapping annotated regions of the sequences as protein sequences
        @param overlapLenMin: the sequences must have at least this overlap
        @param overlapAnnotLenMin: the overlapping annotated regions must have at least this overlap
        @param stopOverlapMaxMismatch: if this percentage of the positions don't overlap, stop and return 3*[None]
        @param continuousOverlap: count the overlap score for dna sequences at once
        (do not consider annotation separately) in this case the last returned value is the annotation length,
        not score)

        @type rec1: read_rec.ReadRec
        @type rec2: read_rec.ReadRec
        @type overlapIdx: int
        @type tripletMap: dict[str,str]
        @type considerProtSeqCmp: bool
        @type overlapLenMin: int
        @type overlapAnnotLenMin: int
        @type stopOverlapMaxMismatch: float
        @param continuousOverlap: bool

        @rtype: (float | None, float | None, float | None)

        @return: (overlap-probability, expected-length-of-matching-entries, expected-length-of-matching-annot-entries)
    """
    dna1 = rec1.dnaSeq
    dna2 = rec2.dnaSeq

    # get start positions withing the dna sequences
    if overlapIdx >= 0:
        r1 = 0
        r2 = overlapIdx
    else:
        r1 = - overlapIdx
        r2 = 0

    # length of the dna overlap
    overlapLen = min(len(dna1) - r1, len(dna2) - r2)

    # overlap is too short, stop
    if overlapLen < overlapLenMin:
        return tuple(3 * [None])

    # maximum mismatch positions within the overlap allowed
    maxMismatch = int(stopOverlapMaxMismatch * float(overlapLen))

    i = 0
    j1 = r1
    j2 = r2
    mismatchCount = 0
    annLen = 0
    annStart1 = None
    annStart2 = None

    # compare all positions within the overlap, stop if too many mismatches
    while i < overlapLen:

        # within the annotated region
        if (rec1.annotStart <= j1 < rec1.annotStart + rec1.annotLen) \
                and (rec2.annotStart <= j2 < rec2.annotStart + rec2.annotLen):
            if annLen == 0:
                annStart1 = j1
                annStart2 = j2
            annLen += 1
            withinAnnot = True
        else:
            withinAnnot = False

        # overlap mismatch
        if dna1[j1] != dna2[j2]:
            if not considerProtSeqCmp:
                mismatchCount += 1
            elif withinAnnot:
                # do the corresponding triplets encode the same aminoacid
                # get the frame offset, the triplet start-position, the triplet
                offset = (j1 - rec1.annotStart) % 3
                start = j1 - offset
                amino1 = tripletMap.get(dna1[start:start+3])

                offset = (j2 - rec2.annotStart) % 3
                start = j2 - offset
                amino2 = tripletMap.get(dna2[start:start+3])

                # triplets encode different aminoacids
                if not(amino1 is not None and amino1 == amino2):
                    mismatchCount += 1
            else:
                mismatchCount += 1

            # there are too many mismatches
            if mismatchCount >= maxMismatch:
                return tuple(3 * [None])
        i += 1
        j1 += 1
        j2 += 1

    # the annotation overlap is too short
    if annLen < overlapAnnotLenMin:
        return tuple(3 * [None])

    # count the overlap only based on the DNA sequences and the whole overlap at once
    if continuousOverlap and not considerProtSeqCmp:
        overlapProb, overlapScore = rec1.qsArray.getOverlapScore(r1, rec2.qsArray, r2, overlapLen)
        overlapScoreAnnot = annLen
    else:
        # compute probability and scores
        overlapProb = 0.
        overlapScore = 0.

        # before annotation overlap
        if r1 == annStart1:
            assert r2 == annStart2
            bLen = 0
        else:
            bLen = annStart1 - r1
            assert bLen > 0
            p, score = rec1.qsArray.getOverlapScore(r1, rec2.qsArray, r2, bLen)
            w = float(bLen) / float(overlapLen)
            overlapProb += w * p
            overlapScore += score

        # annot overlap
        assert annStart1 is not None and annStart2 is not None and annLen >= 3

        if considerProtSeqCmp:
            # compute the overlap probabilities, within annotated regions consider PROT sequences
            p, overlapScoreAnnot = rec1.qsArray.getOverlapScoreProt(annStart1, rec2.qsArray, annStart2, annLen)
        else:
             # compute the overlap probabilities considering only the dna sequences
            p, overlapScoreAnnot = rec1.qsArray.getOverlapScore(annStart1, rec2.qsArray, annStart2, annLen)

        w = float(annLen) / float(overlapLen)
        overlapProb += w * p
        overlapScore += overlapScoreAnnot

        if bLen + annLen == overlapLen:
            aLen = 0
        else:
            aLen = overlapLen - bLen - annLen
            assert aLen > 0
            p, score = rec1.qsArray.getOverlapScore(annStart1 + annLen, rec2.qsArray, annStart2 + annLen, aLen)
            w = float(aLen) / float(overlapLen)
            overlapProb += w * p
            overlapScore += score

        assert overlapLen == bLen + annLen + aLen

    return (overlapProb, overlapScore, overlapScoreAnnot)


def getPossibleOverlaps(rec1, rec2):
    """
        Get possible overlaps of the read-records.

        @type rec1: read_rec.ReadRec
        @type rec1: read_rec.ReadRec

        @return: a list of possible overlaps as index 0 of rec1 start within rec2 (negative values possible, also None)
        @rtype: list[int] | None
    """
    # get the overlap length of the Hmm annotated regions (startHmmCoord + local annotation), this is just an estimate
    start = max(3 * rec1.hmmCoordStart, 3 * rec2.hmmCoordStart)
    end = min(3 * rec1.hmmCoordStart + rec1.annotLen - 1, 3 * rec2.hmmCoordStart + rec2.annotLen - 1)

    # there is "approximately" an overlap of the Hmm annotated regions
    if start <= end:
        overlapList = []

        # difference of the Hmm coordinates start points
        coordDiff = 3 * (rec1.hmmCoordStart - rec2.hmmCoordStart)

        # pos within rec2 matching start of rec1
        middlePoint = rec2.annotStart + coordDiff - rec1.annotStart

        # how far to go to the left and right from the middle point
        goLeftMax = rec1.annotLen + coordDiff  # this may be reduced to boost the performance!
        goRightMax = rec2.annotLen - coordDiff  # may be reduced as well to boost the performance!

        overlapList.append(middlePoint)

        # go left and right and add all possible overlaps
        diff = 3
        goLeft = True
        goRight = True
        while goLeft or goRight:
            if diff < goLeftMax:
                overlapList.append(middlePoint - diff)
            else:
                goLeft = False
            if diff < goRightMax:
                overlapList.append(middlePoint + diff)
            else:
                goRight = False
            diff += 3

        return overlapList

    else:
        # there is no overlap
        return None


def runSnowball(recList, recSeedList, considerProtSeqCmp=False, considerOnlyPOverlap=False,
                pOverlapMin=0.8, scoreOverlapLenMin=0.8, scoreOverlapAnnotLenMin=0.2,
                stopOverlapMaxMismatch=0.1, translTable=11):
    """
        Runs the snowball assembly algorithm. Given read records representing individual reads,
        outputs super-reads (i.e. contigs).

        @param recList: read-record list sorted according to the start alignment coordinates within a gene family
        @param recSeedList: read-record list sorted according to the biggest overlap with all the annotations

        @param considerProtSeqCmp: annotated regions are considered as PROT for overlap score calculations
        @param considerOnlyPOverlap: join two (super-)reads only based on the probability scores
        @param pOverlapMin: min. overlap probability required
        @param scoreOverlapLenMin: min. overlap expected length
        @param scoreOverlapAnnotLenMin: min. overlap annotation expected length
        @param stopOverlapMaxMismatch: stop investigating an overlap when so many positions DNA have a mismatch
        @param translTable: 11 for bacteria and archaea

        @type recList: list[read_rec.ReadRec]
        @type recSeedList: list[read_rec.ReadRec]
        @type considerProtSeqCmp: bool
        @type considerOnlyPOverlap: bool
        @type pOverlapMin: float
        @type scoreOverlapLenMin: float
        @type scoreOverlapAnnotLenMin: float
        @type stopOverlapMaxMismatch: float
        @type translTable: int

        @return: set of read-records representing assembled super-reads (contigs)
        @rtype: set[read_rec.ReadRec]
    """
    if len(recList) == 0 or len(recSeedList) == 0:
        return set()

    # initialize the working set by the seed
    seed = recSeedList[0]
    workingSet = {seed}

    # map: triplet -> aminoacid
    tripletMap = com.getTripletMap(translTable)

    # next record to be considered in the snowball
    left = None
    right = None

    # find the position of the seed in the recList
    for i in com.binarySearch(recList, seed, fun=lambda x: x.hmmCoordStart):
        if recList[i].recordId == seed.recordId:
            if i - 1 >= 0:
                left = i - 1
            if i + 1 < len(recList):
                right = i + 1
            break

    # while there is a record to be considered for joining with one of the seeds in the working set
    while (left is not None) or (right is not None):

        # list of pairs (seed, left or right) to be inspected for joining
        inspectList = []

        for seed in workingSet:
            if left is not None:
                inspectList.append((seed, recList[left]))
            if right is not None:
                inspectList.append((seed, recList[right]))

        # left and right may be the best to join
        if left is not None and right is not None:
            inspectList.append((recList[left], recList[right]))

        # get the best overlap, compute overlap of the left and right to all the seeds (and left and right)
        joinInfo = getBestOverlap(inspectList, tripletMap, considerProtSeqCmp, pOverlapMin, scoreOverlapLenMin,
                                  scoreOverlapAnnotLenMin, stopOverlapMaxMismatch, considerOnlyPOverlap)

        moveRight = False
        moveLeft = False
        # there is a sufficient overlap of the seed and a neighbor, join them
        if joinInfo is not None:

            # remove the seed from the working set (if the best overlap is not between left and right)
            if not (left is not None and joinInfo.seedRec.recordId == recList[left].recordId):
                workingSet.remove(joinInfo.seedRec)

            # merge the records (seed and a neighbor)
            mergedRec = read_rec.mergeRecords(joinInfo.seedRec, joinInfo.neighborRec, joinInfo.overlap)

            # replace the seed in the working set by the merged record
            workingSet.add(mergedRec)

            # was the neighbor left or right
            if left is not None and joinInfo.neighborRec.recordId == recList[left].recordId:
                moveLeft = True
            else:
                assert right is not None and joinInfo.neighborRec.recordId == recList[right].recordId
                moveRight = True

            # if left and right were joined then seed ~ left record
            if left is not None and joinInfo.seedRec.recordId == recList[left].recordId:
                moveLeft = True

        else:
            # there is no sufficient overlap, add neighbors to the working set
            if left is not None:
                workingSet.add(recList[left])
                moveLeft = True
            if right is not None:
                workingSet.add(recList[right])
                moveRight = True

        # shift the indices to get new neighbors to consider
        if moveRight:
            if right + 1 < len(recList):
                right += 1
            else:
                right = None
        if moveLeft:
            if left - 1 >= 0:
                left -= 1
            else:
                left = None

    return workingSet
