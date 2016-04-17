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

    I/O handler. Module version 1.2
"""

import os
import gzip
import cPickle

from algbioi.com import fq
from algbioi.com import rand
from algbioi.com import fasta as fas
from algbioi.haplo import qs as qs_man
from algbioi.haplo import read_rec
from algbioi.haplo import hmm
from algbioi.hsim import comh


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
    nameToDom = hmm.readDomblout(inDomtblout)

    assert inProtFna is None or len(nameToProtSeq) == len(nameToDom)

    # read in pair-end reads, create ReadRec
    for readName, dna, comment, qs in fq.ReadFqGen(inFq):

        readName = readName[1:]  # strip starting @

        protSeq = nameToProtSeq.get(readName, None)

        hit, frameTag = nameToDom[readName]

        annotStart, annotLen, strain, score, acc = hmm.dnaHitInfo(hit, dna, protSeq)  # strain, score, acc not used

        if strain == 1:
            assert 1 <= frameTag <= 3
        else:
            assert 4 <= frameTag <= 6

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

        tokens = comment.split('\t')

        if len(tokens) == 5:
            # get the ends of the pair-end read and corresponding quality-scores
            p, dna1, qs1, dna2, qs2 = tokens

            # get the QSArray representation of the first read-end
            qsA1 = qs_man.QSArray(dna=dna1, qsArrayFq=qs1)

            # get the reverse complement of the second read-end
            qsA2 = qs_man.QSArray(dna=dna2, qsArrayFq=qs2)
            qsA2.revCompl()

            # get the consensus QS Array representing of the joined read
            qsA = qs_man.QSArray(qsA1=qsA1, qsA2=qsA2, pos1Within2=len(dna2) - len(dna))
        else:
            # there is just a simple dna sequence (i.e. not joined reads)
            qsA = qs_man.QSArray(dna=dna, qsArrayFq=qs)

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


def getReadNameToHitList(domFile):
    """
        map: read-name -> list of hits (lists sorted according to the scores)
    """
    nameToHitList = {}
    for line in gzip.open(domFile):
        line = line.strip()
        if line.startswith('#'):
            continue
        assert line.startswith('@')
        tokens = line.split()
        name = tokens[0][1:-2]

        if name not in nameToHitList:
            nameToHitList[name] = []
        nameToHitList[name].append(tokens)

    # sort lists according to the best hit score
    for hitList in nameToHitList.values():
        hitList.sort(key=lambda x: float(x[13]), reverse=True)

    return nameToHitList


def partitionReads(sampleDir, scoreThreshold, accuracyThreshold, shuffleRandSeed, pfamPartitionedDir, joinedReads=True,
                   considerSam=True):
    """
        Partitioning reads into the individual gene-domains.
    """
    try:
        strainDirList = map(lambda x: os.path.join(sampleDir, x), os.listdir(sampleDir))

        samplePartDir = os.path.join(sampleDir, pfamPartitionedDir)
        if not os.path.isdir(samplePartDir):
            os.mkdir(samplePartDir)

        # for each gene-dom
        fqOutDict = {}
        fqProtOutDict = {}
        fqSamOutDict = {}
        fqDomOutDict = {}

        for strainDir in strainDirList:
            strainAcc = os.path.basename(strainDir)
            if strainAcc == 'sample_partitioned':  # skip the directory containing the partitioned data
                continue
            if os.path.isdir(strainDir):
                # for each dom file
                for f in os.listdir(strainDir):

                    if (joinedReads and f.endswith('join_prot.domtblout.gz')) \
                            or (not joinedReads and f.endswith('pair1_prot.domtblout.gz')):  # a domtblout file found
                        i = f.split('_', 1)[0]
                        if i.isdigit():
                            if joinedReads:
                                domPath = os.path.join(strainDir, f)
                                domPath2 = None
                                fqPath = os.path.join(strainDir, '%s_join.fq.gz' % (i,))
                                fqPair1Path = os.path.join(strainDir, '%s_pair1.fq.gz' % (i,))
                                fqPair2Path = os.path.join(strainDir, '%s_pair2.fq.gz' % (i,))
                                fqProtPath = os.path.join(strainDir, '%s_join_prot.fna.gz' % (i,))
                                fqProtPath2 = None
                                samPath = os.path.join(strainDir, '%s_join_gmap.sam.gz' % (i,))
                                assert os.path.isfile(fqPath)
                            else:
                                domPath = os.path.join(strainDir, f)
                                domPath2 = os.path.join(strainDir, '%s_pair2_prot.domtblout.gz' % (i,))
                                fqPath = None
                                fqPair1Path = os.path.join(strainDir, '%s_pair1.fq.gz' % (i,))
                                fqPair2Path = os.path.join(strainDir, '%s_pair2.fq.gz' % (i,))
                                fqProtPath = os.path.join(strainDir, '%s_pair1_prot.fna.gz' % (i,))
                                fqProtPath2 = os.path.join(strainDir, '%s_pair2_prot.fna.gz' % (i,))
                                samPath = os.path.join(strainDir, '%s_pair.sam.gz' % (i,))
                                assert os.path.isfile(domPath2) and os.path.isfile(fqProtPath2)
                            assert os.path.isfile(domPath) and os.path.isfile(fqPair1Path) \
                                   and os.path.isfile(fqPair2Path) and os.path.isfile(fqProtPath)

                            if considerSam:
                                assert os.path.isfile(samPath)

                            # map: read-name -> list of hits (lists sorted according to the scores)
                            nameToHitList = getReadNameToHitList(domPath)
                            if not joinedReads:
                                nameToHitList2 = getReadNameToHitList(domPath2)
                                len1 = len(nameToHitList)
                                len2 = len(nameToHitList2)
                                nameToHitList.update(nameToHitList2)
                                assert len(nameToHitList) == len1 + len2

                            # map: read-name-prot -> seq-prot
                            protNameToSeq = fas.fastaFileToDictWholeNames(fqProtPath)
                            if not joinedReads:
                                protNameToSeq.update(fas.fastaFileToDictWholeNames(fqProtPath2))

                            # map: read-name -> sam-line-entry
                            if considerSam:
                                readNameToSam = {}
                                if joinedReads:
                                    for line in gzip.open(samPath):
                                        line = line.strip()
                                        if line.startswith('#'):
                                            continue
                                        readName = line.split('\t', 1)[0]
                                        # lines with only 11 entries will be padded with * to 12
                                        if len(line.split('\t')) == 11:
                                            line += '\t*'
                                        readNameToSam[readName] = line + '\t' + strainAcc
                                else:
                                    entry = []
                                    for line in gzip.open(samPath):
                                        line = line.strip()
                                        if line.startswith('#') or line.startswith('@'):
                                            continue
                                        if len(entry) < 2:
                                            entry.append(line)
                                        if len(entry) == 2:
                                            readName = entry[0].split('\t', 1)[0]
                                            assert readName == entry[1].split('\t', 1)[0]
                                            readNameToSam[readName] = entry[0] + '\t*\t' + strainAcc + '\n' \
                                                                      + entry[1] + '\t*\t' + strainAcc
                                            entry = []

                            # map read-name -> "pair1-dna tab pair1-qs tab pair2-dna tab pair2-qs"
                            if joinedReads:
                                readNameToPairReads = fq.getFqToDict(fqPair1Path, fqPair2Path)
                            else:
                                readNameToPairReads = None

                            if joinedReads:
                                g1 = fq.ReadFqGen(fqPath)
                                g2 = []
                            else:
                                g1 = fq.ReadFqGen(fqPair1Path)
                                g2 = fq.ReadFqGen(fqPair2Path)

                            # go over all reads
                            for readName, dna, p, qs in list(g1) + list(g2):
                                readName = readName[1:]  # strip starting '@'

                                # take the hit with the highest score
                                topHit = None
                                if readName in nameToHitList:
                                    topHit = nameToHitList[readName][0]

                                # is the hit significant, filter according to the score and accuracy
                                if topHit is None or float(topHit[13]) < scoreThreshold or float(topHit[21]) < accuracyThreshold:
                                    continue
                                else:
                                    famName = topHit[3]
                                    if famName not in fqOutDict:
                                        fqOutDict[famName] = []
                                        fqProtOutDict[famName] = []
                                        fqSamOutDict[famName] = []
                                        fqDomOutDict[famName] = []

                                    if joinedReads:
                                        comment = readNameToPairReads[readName]
                                    else:
                                        comment = ''
                                    fqOutDict[famName].append((readName, dna, qs, comment))

                                    protSeqName = topHit[0]
                                    protSeq = protNameToSeq[protSeqName]
                                    fqProtOutDict[famName].append((readName, protSeq))

                                    if considerSam:
                                        if joinedReads:
                                            fqSamOutDict[famName].append(readNameToSam[readName])
                                        else:
                                            fqSamOutDict[famName].append(readNameToSam[readName[:-2]])

                                    # top hit coordinates within the read
                                    startOnRead, overlapLen, strain = hmm.dnaHitInfo(topHit, dna, protSeq)[:3]

                                    fqDomOutDict[famName].append('\t'.join(topHit) + '\t%s\t%s\t%s' % (startOnRead, overlapLen, strain))

        if joinedReads:
            ident = 'join'
        else:
            ident = 'pair'

        # for each gene dom, store reads into a file
        for famName, fqContentList in fqOutDict.iteritems():

            # get the tagged fam-dom-name that can be used in file names
            pf = comh.getGeneNameToFileName(famName)[:-4]

            # define output files
            fqOutO = os.path.join(samplePartDir, 'o_%s_%s.fq.gz' % (pf, ident))  # 'o_' ~ ordered entries
            fqOutR = os.path.join(samplePartDir, 'r_%s_%s.fq.gz' % (pf, ident))  # 'r_' ~ random shuffled entries

            fqProtOutO = os.path.join(samplePartDir, 'o_%s_%s_prot.fna.gz' % (pf, ident))
            fqProtOutR = os.path.join(samplePartDir, 'r_%s_%s_prot.fna.gz' % (pf, ident))

            fqSamOutO = os.path.join(samplePartDir, 'o_%s_%s_gmap.sam.gz' % (pf, ident))
            fqSamOutR = os.path.join(samplePartDir, 'r_%s_%s_gmap.sam.gz' % (pf, ident))

            fqDomOutO = os.path.join(samplePartDir, 'o_%s_%s_prot.domtblout.gz' % (pf, ident))
            fqDomOutR = os.path.join(samplePartDir, 'r_%s_%s_prot.domtblout.gz' % (pf, ident))

            # write FASTQ
            fqOut = fq.WriteFq(fqOutO)
            for e in fqContentList:
                fqOut.writeFqEntry('@' + e[0], e[1], e[2], e[3])
            fqOut.close()

            # write PROT
            fqProtOut = fq.WriteFq(fqProtOutO)
            fqProtOut.write('\n'.join(map(lambda x: '>%s\n%s' % (x[0], x[1]), fqProtOutDict[famName])) + '\n')
            fqProtOut.close()

            # write SAM
            if considerSam:
                fqSamOut = fq.WriteFq(fqSamOutO)
                fqSamOut.write('\n'.join(fqSamOutDict[famName]) + '\n')
                fqSamOut.close()

            # write DOM
            fqDomOut = fq.WriteFq(fqDomOutO)
            fqDomOut.write('\n'.join(fqDomOutDict[famName]) + '\n')
            fqDomOut.close()

            # shuffle file entries (to remove any bias imposed by the ordering)
            rand.shuffleLines(fqOutO, fqOutR, 4, shuffleRandSeed)
            rand.shuffleLines(fqProtOutO, fqProtOutR, 2, shuffleRandSeed)
            rand.shuffleLines(fqDomOutO, fqDomOutR, 1, shuffleRandSeed)

            if considerSam:
                if joinedReads:
                    rand.shuffleLines(fqSamOutO, fqSamOutR, 1, shuffleRandSeed)
                else:
                    rand.shuffleLines(fqSamOutO, fqSamOutR, 2, shuffleRandSeed)

            # delete ordered files (keep only the shuffled ones)
            os.remove(fqOutO)
            os.remove(fqProtOutO)
            if considerSam:
                os.remove(fqSamOutO)
            os.remove(fqDomOutO)
    except Exception as e:
        print('Exception in partitionReads:')
        print sampleDir, scoreThreshold, accuracyThreshold, shuffleRandSeed, pfamPartitionedDir, joinedReads
        print e.message
        print type(e)
        print e.args
        raise e

# TESTS ---------------------------------------------------

def _test():
    baseDir = '/home/igregor/Documents/work/hsim/562/samples/0/sample_partitioned'
    baseName = 'aminotran_3_1'  # 'duf2158_3'
    recList = parse(inFq=os.path.join(baseDir, 'r_' + baseName + '_join.fq.gz'),
                    inDomtblout=os.path.join(baseDir, 'r_' + baseName + '_join_prot.domtblout.gz'),
                    inProtFna=os.path.join(baseDir, 'r_' + baseName + '_join_prot.fna.gz'))
    f = '/home/igregor/Documents/work/hsim/562/samples/0/tmp_test.gz'
    storeReadRec(recList, f)
    recList2 = loadReadRec(f)
    for r1, r2 in zip(recList, recList2):
        assert r1.dnaSeq == r2.dnaSeq
        for i in range(len(r1.posCovArray)):
            assert r1.posCovArray[i] == r2.posCovArray[i]
    print len(recList)
    print len(recList2)
