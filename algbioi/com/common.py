"""
    Common functionality module version: 1.2

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
"""

import os
import re
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def getTripletMap(translTable=11):
    """
        Gets a mapping between nucleotide triplets and aminoacids

        @return: map: dna-triplet -> aminoacid
        @rtype: dict[str, str]
    """
    nucL = ['A', 'T', 'G', 'C']
    tripletToAmino = {}
    for p1 in nucL:
        for p2 in nucL:
            for p3 in nucL:
                dna = p1 + p2 + p3
                prot = str(Seq(dna, generic_dna).translate(table=translTable))
                assert dna not in tripletToAmino
                tripletToAmino[dna] = prot

    return tripletToAmino


def isclose(a, b, epsilon=0.00000001):
    """
        Analogy to numpy.isclose(a, b)
    """
    try:
        return np.isclose(a, b)
    except AttributeError:
        if abs(a - b) < epsilon:
            return True
        else:
            return False


def removeNonDna(seq):
    """
        Replaces stretches of nonDNA characters with one 'N'.
    """
    return re.sub(r'[^ATGCatgc]+', 'N', seq).upper()


def noNewLine(line):
    """
        Delete all '\n' and '\r' characters in a string.
    """
    return line.replace('\n', '').replace('\r', '')


def createTagFilePath(dstDir, fileNameFromPath, tag):
    """
        Returns a path that results from the concatenation of dstDir and a file from filePath where the name of the file
        changes from e.g.: input.n.fas" to "input.n._ids.fas"
    """
    lastDotIdx = fileNameFromPath.rfind('.')
    return os.path.join(os.path.normpath(dstDir),
                        os.path.basename(os.path.normpath(str(fileNameFromPath[0:lastDotIdx]
                                                              + fileNameFromPath[lastDotIdx:]
                                                              + '.' + tag))))


def getMothurOutputFilePath(inputFastaFilePath, refTaxonomyFilePath, suffix='.taxonomy'):
    """
        Returns mothur prediction file path that is generated based on the arguments of the Mothur classify command.
    """
    dirName = os.path.dirname(inputFastaFilePath)
    fastaBaseName = os.path.basename(inputFastaFilePath)
    fastaPart = fastaBaseName[0:fastaBaseName.rindex('.')]  # without suffix
    taxBaseName = os.path.basename(refTaxonomyFilePath)
    taxPart = taxBaseName[0:taxBaseName.rindex('.')]  # without suffix
    if '.' in taxPart:
        taxPart = taxPart[(taxPart.rindex('.') + 1):]  # from last comma till the end

    return os.path.join(dirName, str(fastaPart + '.' + taxPart + suffix))


def seqFileCmp(file1, file2):
    """
        Returns true if two files contain the same sequences regardless of their names (and empty spaces).
    """
    seqList1 = seqFileToSeqList(file1)
    seqList2 = seqFileToSeqList(file2)
    if len(seqList1) != len(seqList2):
        print "The files contain different number of sequences", len(seqList1), len(seqList2)
        return False
    else:
        print "Number of sequences: ", len(seqList1)

    seqFoundInS2 = set()
    seqAF = set()
    idx1 = 0
    for s1 in seqList1:
        idx1 += 1
        idx2 = 0
        for s2 in seqList2:
            idx2 += 1
            if s1 == s2:
                #print "same:", idx1, idx2
                if s2 not in seqAF:
                    #print idx1, idx2
                #else:
                    seqAF.add(s2)
                if idx2 not in seqFoundInS2:
                    #print "One sequence is in one file more than once", idx1, idx2
                    seqFoundInS2.add(idx2)
                    continue

    s1 = set()
    s2 = set()
    for s in seqList1:
        s1.add(s)
    for s in seqList2:
        s2.add(s)
    if len(s1) != len(s2):
        print "The length of unique sequences differ S1:", len(s1), "S2:", len(s2)
    else:
        print "The number of unique sequences is: ", len(s1)

    if len(seqFoundInS2) == len(seqList1):
        print "Both files contain the same sequences"
        return True
    else:
        print "Sequences matches: ", len(seqList1), " Sequences found: ", len(seqFoundInS2)
        return False


def seqFileToSeqList(fileName):
    """
        @deprecated: use functionality of algbioi.com.fasta
    """
    seqList = []
    f = None
    try:
        f = open(os.path.normpath(fileName), 'r')
    except Exception:
        print "Cannot open file:", fileName
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqList.append(seq)  # store seq
                    seq = ''
                name = line.replace('>', '')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqList.append(seq)  # store seq
        return seqList
    finally:
        if f is not None:
            f.close()


class NodeNewick():
    def __init__(self, label, nodeList=None):
        """
            @attention: needs to be tested !!!
            @type nodeList: list of NodeNewick
            @type label: str
        """
        assert nodeList is None or len(nodeList) > 0
        self.label = label
        self.nodeList = nodeList

    def isLeaf(self):
        if self.nodeList is None:
            return True
        return False

    def getChildList(self):
        assert self.nodeList is not None
        return self.nodeList


def getNewick(node):
    """
        Get a tree in a newick format, call with a rood of the tree.

        @attention: needs to be tested !!!
        @type node: NodeNewick
    """
    if node.isLeaf():
        return node.label
    else:
        childNewickList = []
        for child in node.getChildList():
            childNewickList.append(getNewick(child))
        return '(' + ','.join(childNewickList) + ')' + node.label


def binarySearch(objList, obj, fun=lambda x: x):
    """
        Implements a simple binary search, returns a list of matching indices.

        @param objList: a list of items
        @param obj: an item we search
        @param fun: a function that extracts a value from an object for the comparison
        @return: a list of indices or an empty list
        @rtype: list[int]
    """
    first = 0
    last = len(objList) - 1
    objFun = fun(obj)

    while first <= last:
        midpoint = (first + last) / 2
        midpointFun = fun(objList[midpoint])

        if midpointFun == objFun:
            while midpoint > 0:
                if fun(objList[midpoint - 1]) == objFun:
                    midpoint -= 1
                else:
                    break
            retIdxList = []
            while midpoint < len(objList):
                if fun(objList[midpoint]) == objFun:
                    retIdxList.append(midpoint)
                    midpoint += 1
                else:
                    break
            assert len(retIdxList) > 0
            return retIdxList
        else:
            if objFun < midpointFun:
                last = midpoint - 1
            else:
                first = midpoint + 1
    return []
