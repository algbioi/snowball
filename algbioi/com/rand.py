"""
    Random Python module version: 1.2

    The MIT License (MIT)

    Copyright (c) 2015  Ivan Gregor

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

    This module contain functions related to randomness.
"""

import gzip
import cPickle
import tempfile
import time
import numpy as np
from algbioi.com import fq


def shuffleLines(inFilePath, outFilePath, linesPerRecord=4, randomSeed=1, recordLimit=1000000):
    """
        Randomly shuffle lines from the input file and output shuffled lines to the output file.

        Algorithm: split the input files into several smaller files, shuffle individual files, randomly merge the files.

        Estimated throughput: 16-18 MB/s

        @attention: each line must end with the EOL character including the last one !!!

        @param inFilePath: input file path
        @param outFilePath: output file path
        @param linesPerRecord: how many lines should be considered as one record (e.g. 4 for FASTQ files)
        @param randomSeed: set the random seed
        @param recordLimit: the maximum number of lines / records that will be shuffled in the main memory at once
    """
    # get an instance of the random number generator
    rand = np.random.RandomState(randomSeed)
    # a list of temporary files
    tmpFileList = []
    # a list of records to be stored in one temporary file
    recordList = []
    # to store one record
    record = []
    # read records one by one
    if inFilePath.strip().endswith('.gz'):
        inFileOpen = gzip.open(inFilePath, mode='r')
    else:
        inFileOpen = open(inFilePath)
    for line in inFileOpen:
        record.append(line)
        if len(record) == linesPerRecord:
            # store one record
            recordList.append(''.join(record))
            record = []
            # if more than limit records read, store its permutation to a temporary file
            if len(recordList) == recordLimit:
                tmpFileList.append(_storePermutationToTmpFile(recordList, rand))
                recordList = []
    # store the rest of the records as a permutation to a temporary file
    if len(recordList) > 0:
        tmpFileList.append(_storePermutationToTmpFile(recordList, rand))

    # open a file for writing
    outFile = fq.WriteFq(outFilePath)
    # write one record to the output from a randomly chosen file in each iteration
    if len(tmpFileList) > 0:
        while True:
            i = rand.randint(0, len(tmpFileList))
            try:
                record = cPickle.load(tmpFileList[i])
                outFile.write(record)
            except EOFError:
                tmpFileList[i].close()
                del tmpFileList[i]
            if len(tmpFileList) == 1:
                while True:
                    try:
                        record = cPickle.load(tmpFileList[0])
                        outFile.write(record)
                    except EOFError:
                        tmpFileList[0].close()
                        break
                break

    outFile.close()


def _storePermutationToTmpFile(recordList, rand):
    """
        Stores a permutation of the record list to a temporary file and return it.
        The file is ready for reading!

        @param recordList: a list of records
        @param rand: instance of the random generator

        @return: a temporary file containing pickled records (one record one pickled entry) (ready for reading)
    """
    # get a temporary file
    tmp = tempfile.TemporaryFile(mode='rb+')
    # randomly permute the list
    rand.shuffle(recordList)
    # write the permuted records to the temporary file
    for item in recordList:
        cPickle.dump(item, tmp, cPickle.HIGHEST_PROTOCOL)
    tmp.flush()
    tmp.seek(0)
    return tmp


def getRandLognormNumbers(listOfCounts, mean=1, sd=2, minVal=1., maxVal=50., randomSeed=1, rand=None):
    """
        Get a list of lists for each count in the input list.

        @param listOfCounts: each count is a number of strains in one sample
        @param mean: mean of the log-norm distribution
        @param sd: standard deviation of the log-norm distribution
        @param minVal: minimum value generated
        @param maxVal: maximum value generated
        @param randomSeed: random seed used for the random number generator (used only if param rand is None !!!)
        @param rand: use this random state numpy.random.RandomState (default None)

        @return: list of lists
    """
    # get an instance of the random number generator
    if rand is None:
        rand = np.random.RandomState(randomSeed)
    retList = []
    # for each count create a list of random log-norm numbers
    for count in listOfCounts:
        entryList = []
        # generate a random list
        while len(entryList) < count:
            for r in rand.lognormal(mean, sd, count - len(entryList)):
                if minVal <= r <= maxVal:
                    entryList.append(r)

        retList.append(entryList)
    return retList


def getRandSubsets(itemList, subsetSizes, maxSubsetsPerSize, randomSeed=1, rand=None, coverFromRests=0.3):
    """
        For a list of items create a list of lists representing random subsets of the input list.

        For each size, random subsets are chosen greedily.

        @param itemList: list of items
        @param subsetSizes: a tuple of subset sizes
        @param maxCountPerSize: maximum number of subsets generated for one subset size
        @param randomSeed: random seed used for the random number generator (used only if param rand is None !!!)
        @param rand: use this random state numpy.random.RandomState (default None)
        @return: list of lists representing subsets
    """
    if rand is None:
        rand = np.random.RandomState(randomSeed)
    sampleList = []
    itemListLen = len(itemList)
    # for each sample size
    for size in subsetSizes:
        if itemListLen >= size:
            # list of indices for sampling
            idxList = range(itemListLen)
            # get samples for one size
            for tmp in range(maxSubsetsPerSize):
                if len(idxList) >= size:
                    # randomly choose indices for a new sample
                    sampleIdx = rand.choice(idxList, size, replace=False)
                    # create a sample
                    sample = []
                    for i in sampleIdx:
                        sample.append(itemList[i])
                        idxList.remove(i)
                    sampleList.append(sample)

                elif len(idxList) > size * coverFromRests:
                    # create a sample from the remaining indices
                    sample = []
                    rand.shuffle(idxList)
                    idxList2 = range(itemListLen)
                    for i in idxList:
                        sample.append(itemList[i])
                        idxList2.remove(i)
                    sampleIdx = rand.choice(idxList2, size - len(sample), replace=False)
                    for i in sampleIdx:
                        sample.append(itemList[i])
                    sampleList.append(sample)
                    break
                else:
                    break
    return sampleList


def strToRandInt(s, upperBound=31727):
    """
        Get a pseudo random integer given a string.
        Each call with the same string result in the same integer.
        @param upperBound:  a big prime number (e.g. 31727)

        @return: integer from [0, upperBound - 1]
    """
    rand = np.random.RandomState(np.array(map(lambda x: ord(x), list(s))))
    return rand.randint(upperBound)

# TESTS ---------------------------------------------------

def _testShuffleLines():
    t = time.time()
    # inF = '/home/user19/test_reads.fq'
    # inF = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AMUV00000000/pair1.fq'
    # inF = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AOEC00000000/pair1.fq'
    inF = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NC_013008/pair1.fq'
    outF = '/Users/ivan/Documents/nobackup/hsim01/562/samples/0/NZ_AMUV00000000/pair1_rand_.fq'
    # outF = '/home/user19/test_reads_shuffle.fq'
    shuffleLines(inF, outF, randomSeed=1, recordLimit=1000000)
    print('Runtime shuffle: %s' % (time.time() - t))
    # t = time.time()
    # fas.cmpSeqFiles(inF, outF)
    # print('Runtime test: %s' % (time.time() - t))


def _test1():
    print getRandLognormNumbers([3, 5, 7])
    print getRandSubsets(['A', 'B', 'C', 'D', 'E'], subsetSizes=(3,4,5), maxSubsetsPerSize=4)
    print strToRandInt('0NC_0130082.05528149')
    print strToRandInt('0NC_0130082.05528149')
    print strToRandInt('0NC_0130083.05528149')
    print strToRandInt('0NC_0130083.05528149')
    print getRandLognormNumbers([500], mean=1, sd=2, minVal=1., maxVal=50., randomSeed=1, rand=None)
