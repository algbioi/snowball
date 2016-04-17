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

    Implements the main function of the package. Module version 1.2
"""

import numpy as np

from algbioi.haplo import hio
from algbioi.haplo.snowball import runSnowball


def getHmmCovArray(recList):
    """
        Returns an array representing the coverage, i.e. how the HMM "dna alignment" is covered by the HMM annot. reads.

        @type recList list[algbioi.haplo.read_rec.ReadRec]

        @rtype: ndarray
    """
    # find out the alignment length considering all read-records
    aliLen = 0
    for rec in recList:
        aliLen = max(aliLen, 3 * rec.hmmCoordStart + rec.annotLen)

    # alignment coverage array
    aliCovArray = np.zeros(aliLen, np.float64)

    # go over all records
    for rec in recList:
        hmmCoordStart = rec.hmmCoordStart * 3
        posCovArray = rec.getPosCovArray()

        # go over all alignment positions, count the coverage
        j = 0
        for i in range(rec.annotStart, rec.annotStart + rec.annotLen):
            aliCovArray[hmmCoordStart + j] += posCovArray[i]
            j += 1

    return aliCovArray


def findSeed(recList, aliCovArray):
    """
        Finds the seed, i.e. a record that covers the most alignment positions with the highest coverage.

        For all annotated positions of the read, sums up the alignment coverage at the corresponding positions.
        Returns a list sorted according this sum (the biggest first).

        @type recList: list[algbioi.haplo.read_rec.ReadRec]
        @type aliCovArray: ndarray

        @return: seed list, seeds with the highest support first
        @rtype: list[algbioi.haplo.read_rec.ReadRec]
    """
    recSumList = []
    # compute the overlap sum of each read-record
    for rec in recList:
        hmmCoordStart = rec.hmmCoordStart * 3
        covSum = 0.
        j = 0
        for i in range(rec.annotStart, rec.annotStart + rec.annotLen):
            covSum += aliCovArray[hmmCoordStart + j]
            j += 1

        assert covSum > 0.
        recSumList.append((rec, covSum))

    # sort according to the overlap sum
    recSumList.sort(key=lambda x: x[1], reverse=True)

    # return only the records
    retList = []
    for rec in recSumList:
        retList.append(rec[0])
    return retList


def buildSuperReads(inFq, inDomtblout, inProtFna=None, outFile=None,
                    considerProtSeqCmp=(False, False), considerOnlyPOverlap=(False, True),
                    pOverlapMin=(0.8, 0.8), scoreOverlapLenMin=(0.5, 0.4),
                    scoreOverlapAnnotLenMin=(0.2, 0.2), stopOverlapMaxMismatch=(0.1, 0.05), maxLoops=1, translTable=11):
    """
        Main function of the snowball algorithm.

        Recommendation:
            pOverlap: 0.7 - 0.8
            score: 0.25 - 0.8
            use default values
    """
    try:
        assert len(considerProtSeqCmp) == len(considerOnlyPOverlap) == len(pOverlapMin) == len(scoreOverlapLenMin) \
               == len(scoreOverlapAnnotLenMin) == len(stopOverlapMaxMismatch) >= maxLoops

        # 1. read in read records
        recList = hio.parse(inFq, inDomtblout, inProtFna)

        # 2. sort read record list according to the HMM start positions
        recList.sort(key=lambda x: x.hmmCoordStart)

        # 3. build the alignment coverage array
        aliCovArray = getHmmCovArray(recList)

        # 4. find the hotspot (starting seed)
        seedList = findSeed(recList, aliCovArray)

        # 5. run the snowball algorithm
        recSet = runSnowball(recList, seedList, considerProtSeqCmp[0], considerOnlyPOverlap[0], pOverlapMin[0],
                             scoreOverlapLenMin[0], scoreOverlapAnnotLenMin[0], stopOverlapMaxMismatch[0], translTable)

        # 6. optionally, the snowball algorithm can iterate,
        # however, we haven't seen significant improvement when doing so
        i = 1
        inLen = len(recSet)
        while i < maxLoops:
            # get the record list again
            recList = list(recSet)

            # sort: longest seq. first
            recList.sort(key=lambda x: len(x.dnaSeq), reverse=True)
            seedList = recList[:]

            # sort according to the HMM start positions
            recList.sort(key=lambda x: x.hmmCoordStart)

            # run Snowball again
            recSet = runSnowball(recList, seedList, considerProtSeqCmp[i], considerOnlyPOverlap[i], pOverlapMin[i],
                                 scoreOverlapLenMin[i], scoreOverlapAnnotLenMin[i], stopOverlapMaxMismatch[i],
                                 translTable)
            outLen = len(recSet)
            if inLen == outLen:
                break
            inLen = outLen
            i += 1

        # Here, a post-processing step could be performed, i.e. trim low coverage (low quality) ends of the contigs.
        # Nevertheless, we decided not to perform such a step here and give the user freedom to perform this step
        # if needed. (I.e. a user can set custom thresholds)

        # 7. store the results
        if outFile is not None:
            hio.storeReadRec(list(recSet), outFile)
            return None
        else:
            return recSet
    except Exception as e:
        print('Exception in buildSuperReads:')
        print inFq, inDomtblout, inProtFna, outFile
        print e.message
        print type(e)
        print e.args
