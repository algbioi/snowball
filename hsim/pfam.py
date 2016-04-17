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


    ***********************************************************************

    Manage Pfam-A, get mapping between bacterial genes and Pfam gene families.
"""

import os
# import sys
import gzip
# import itertools
# import operator as op
import numpy as np
from algbioi.com import parallel
from algbioi.com import fasta as fas
from algbioi.com import csv
from algbioi.com import fq
from algbioi.com import rand
from algbioi.com import common as com

from algbioi.hsim import comh

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def _main():

    for spec in comh.SPECIES_LIST:
        specDir = os.path.join(comh.REFERENCE_DIR_ROOT, spec)

        # translate gene DNA to gene PROT sequences
        if False:
            dnaToProt(specDir)

        # run HMM search to match gene PROT sequences to PFAM HMMs
        if False:
            searchForGeneFam(specDir)

        # get mapping between geneName and domainName
        if False:
            getMapGeneNameAndPfamDomainName(specDir)

        # get ambiguous domains
        if False:
            getAmbiguousDomains(specDir)

    # partition the file containing Hmm profiles into individual files, one profile per file
    if False:
        getPartitionHmmProfiles(os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILE),
                                os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILES_PARTITIONED_DIR))


def dnaToProt(specDir):
    """
        Translate gene DNA sequences to PROT gene sequences, create one file per geneName.
        Uses translation table 11.
    """
    print('Translating gene DNA sequences to PROT sequences.')

    srcGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_DIR_NAME)
    dstProtGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)

    # create dst directory
    if not os.path.isdir(dstProtGeneDir):
        os.mkdir(dstProtGeneDir)

    # translate sequences
    for f in os.listdir(srcGeneDir):
        fas.dnaToProt(os.path.join(srcGeneDir, f), os.path.join(dstProtGeneDir, f), translTable=11)


def searchForGeneFam(specDir):
    """
        Running hmmsearch to match genes and Pfam domains.
        Very time consuming!
    """
    print('Searching for Pfam gene families and AMPHORA2 genes')
    srcProtGeneDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)
    taskList = []
    # for each gene DNA fasta file
    for f in os.listdir(srcProtGeneDir):
        filePath = os.path.join(srcProtGeneDir, f)
        basePath = filePath
        if filePath.endswith('.fna'):
            basePath = filePath[:-4]
        # define hmmsearch command
        cmd = '%s/hmmsearch -o /dev/null --noali --domtblout %s -E %s --cpu %s %s %s' % \
              (comh.HMMER_BINARY,
               # str(pathBase + '.out'),
               # str(pathBase + '.aln'),  # -A %s
               # str(pathBase + '.tblout'),  # --tblout %s
               str(basePath + '.domtblout'),
               0.01,
               1,
               os.path.join(comh.HMM_PROFILE_DIR, comh.HMM_PROFILE_FILE), filePath)
        taskList.append(parallel.TaskCmd(cmd, cwd=srcProtGeneDir))

    # run commands in parallel
    parallel.reportFailedCmd(parallel.runCmdParallel(taskList, comh.MAX_PROC))


def getMapGeneNameAndPfamDomainName(specDir):
    """
        Map Pfam-domain names onto geneNames.
    """
    print("Mapping Pfam-domain names onto gene names.")
    srcDir = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME)
    assert os.path.isdir(srcDir)

    # output mapping files
    geneNameToDomOut = csv.OutFileBuffer(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_GENE_NAME_TO_PFAM))
    domToGeneNameOut = csv.OutFileBuffer(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME))

    geneNameToDomEntryList = {}

    # for each dom file
    for f in os.listdir(srcDir):
        if f.endswith('.domtblout'):
            seqNameToList = {}
            geneName = None

            # for each line in the dom file
            for line in open(os.path.join(srcDir, f)):
                # skip comments
                if line.startswith('#'):
                    continue

                # parse a line
                tokens = line.split()
                assert len(tokens) == 23, '%s, %s, %s' % (f, tokens, len(tokens))

                targetName, accession1, tlen, queryName, accession2, qlen, eValue, score1, bias, r1, r2, cEvalue, \
                    iEvalue, score2, bias, from1, to1, from2, to2, from3, to3, acc, descriptionOfTarget = tokens

                # get the gene name
                if geneName is None:
                    headerD = dict(zip(map(lambda x: x.split(':')[0], targetName.split(';')),
                                       map(lambda x: x.split(':')[1], targetName.split(';'))))
                    geneName = headerD['geneName']

                # store the entry: gene name ->  (score, dom name)
                entry = (float(score1), queryName)
                if targetName not in seqNameToList:  # first record for the gene sequence name
                    seqNameToList[targetName] = [entry]
                else:
                    seqNameToList[targetName].append(entry)  # there are more records for the gene sequence name

            # for each sequence name, find the match with the highest score
            for targetName in seqNameToList:
                seqNameToList[targetName].sort(key=lambda x: x[0], reverse=True)

            # get a list of dom names
            domNameList = []
            for entry in seqNameToList.values():
                domNameList.append(entry[0])

            # store map: gene-name -> list of (score, dom-name)
            if len(domNameList) > 0:
                geneNameToDomEntryList[geneName] = domNameList

    # map: dom-name -> list of (gene-name, score)
    domNameToGeneNameList = {}

    # distribution: geneName -> count of dom names
    geneNameDist = {}

    # store map: gene-name -> dom-name list
    for geneName, domEntryList in geneNameToDomEntryList.iteritems():
        entrySet = set(map(lambda x: x[1], domEntryList))

        # get the support for each dom-name
        domToSupport = {}
        scoreSum = 0.
        for e in domEntryList:
            domName = e[1]
            score = e[0]
            scoreSum += score
            if domName in domToSupport:
                domToSupport[domName] += score
            else:
                domToSupport[domName] = score
        for domName in domToSupport.keys():
            domToSupport[domName] = (domToSupport[domName] / scoreSum) * 100.
        domNameTupleList = []
        for dom, support in domToSupport.iteritems():
            domNameTupleList.append((dom, support))
        domNameTupleList.sort(key=lambda x: x[1], reverse=True)
        outList = []
        for k, v in domNameTupleList:
            outList.append(k)
            outList.append(str(round(v, 1)))
        # store map: gene-name -> list of (dom-name, support)
        geneNameToDomOut.writeText('%s\t%s\n' % (geneName, ','.join(outList)))
        # update distribution
        dLen = len(entrySet)
        if dLen not in geneNameDist:
            geneNameDist[dLen] = 1
        else:
            geneNameDist[dLen] += 1
        # update mapping: dom name -> gene name
        for scoreDom in domEntryList:
            if scoreDom[1] not in domNameToGeneNameList:
                domNameToGeneNameList[scoreDom[1]] = []
            domNameToGeneNameList[scoreDom[1]].append((geneName, scoreDom[0]))

    # map: dom-name -> list of (gene-name, score)
    domNameToGeneScoreList = {}
    for domName, entryList in domNameToGeneNameList.iteritems():
        geneToScore = {}
        scoreSum = 0.
        for entry in entryList:
            geneName = entry[0]
            score = entry[1]
            scoreSum += score
            if geneName in geneToScore:
                geneToScore[geneName] += score
            else:
                geneToScore[geneName] = score
        geneScoreTupleList = []
        for geneName in geneToScore.keys():
            geneScoreTupleList.append((geneName, (geneToScore[geneName] / scoreSum) * 100.))

        geneScoreTupleList.sort(key=lambda x: x[1], reverse=True)

        domNameToGeneScoreList[domName] = []
        for geneScore in geneScoreTupleList:
            domNameToGeneScoreList[domName].append(geneScore[0])
            domNameToGeneScoreList[domName].append(str(geneScore[1]))

    # distribution: dom name -> count of gene names
    domNameDist = {}

    # store map: dom name -> gene name
    for domName, geneNameScoreList in domNameToGeneScoreList.iteritems():
        domToGeneNameOut.writeText('%s\t%s\n' % (domName, ','.join(geneNameScoreList)))
        # update distribution
        gLen = len(geneNameScoreList) / 2
        if gLen not in domNameDist:
            domNameDist[gLen] = 1
        else:
            domNameDist[gLen] += 1

    geneNameToDomOut.close()
    domToGeneNameOut.close()
    print('geneName -> domName distr: %s' % str(geneNameDist))
    print('domName -> geneName distr: %s' % str(domNameDist))


def getAmbiguousDomains(specDir):
    """
        For Pfam families to which more than one gene-name maps, check, whether both (all) of the gene-names are
        contained in one genome. Add a third column!
    """
    print('Getting ambiguous domains')
    out = csv.OutFileBuffer(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME_AMBIGUOUS))
    totalDomains = 0
    domainAmbiguous = 0
    # for each mapping: dom-name -> list of (gene-name, score)
    for line in open(os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME)):
        line = line.strip()
        domName, geneList = line.split('\t')[:2]
        geneList = geneList.split(',')
        ambiguous = False

        # this dom-name maps to at least two gene-names
        if len(geneList) > 3:
            geneToScore = {}
            for i in range(len(geneList) / 2):
                geneToScore[geneList[2*i]] = float(geneList[2*i+1])

            # map: gene-name -> set of (draft) genomes reference accessions
            geneToRefSet = {}
            geneList = []
            for gene in geneToScore.keys():
                geneList.append(gene)
                filePath = os.path.join(specDir, comh.FASTA_PULL_GENES_PHYLO_DIR_NAME,
                                        comh.getGeneNameToFileName(gene))
                accSet = set()
                for header in fas.fastaFileToDictWholeNames(filePath).keys():
                    accSet.add(comh.getSeqNameToDict(header)['accession'])
                geneToRefSet[gene] = accSet

            # find intersections, two genes from one dom mapped to the same reference
            for i in range(len(geneList) - 1):
                gene1 = geneList[i]
                for gene2 in geneList[i+1:]:
                    # two genes maps to the same reference
                    if not geneToRefSet[gene1].isdisjoint(geneToRefSet[gene2]):
                        # write entry to the file
                        out.writeText('%s\t%s\t%s\t%s\n' % (domName, gene1, gene2,
                                                    ','.join(geneToRefSet[gene1].intersection(geneToRefSet[gene2]))))
                        ambiguous = True
        if ambiguous:
            domainAmbiguous += 1
        totalDomains += 1
    print('Total: %s, Ambiguous: %s' % (totalDomains, domainAmbiguous))
    out.close()


def getPartitionHmmProfiles(srcHmmProfileFile, dstDir):
    """
        Partition a file containing hmm profiles into individual files, each containing one hmm profile.

        @param srcHmmProfileFile: file containing hmm profiles in the stockholm format
        @param dstDir: directory containing individual files
        @type srcHmmProfileFile: str
        @type dstDir: str
    """
    assert os.path.isfile(srcHmmProfileFile)
    assert os.path.isdir(os.path.dirname(dstDir))
    if not os.path.isdir(dstDir):
        os.mkdir(dstDir)
    buff = []
    name = None
    for line in open(srcHmmProfileFile):
        if line.startswith('//'):
            buff.append(line)
            assert name is not None
            out = csv.OutFileBuffer(os.path.join(dstDir, comh.getGeneNameToFileName(name)[:-4]) + '.hmm')
            out.writeText(''.join(buff))
            out.close()
            buff = []
            name = None
        else:
            if line.startswith('NAME'):
                assert name is None
                name = line.strip().split()[1]
            buff.append(line)


class FamInfoHolder():
    def __init__(self, geneFamily):
        self.geneFamily = geneFamily

        # total bp annotated to this family even though there is no gene
        self.ngAnnotatedTotalBp = 0.  # nf.. no family (bp annotated to no family)

        # correct family annotation
        self.cfCa = 0.  # ca.. correctly annotated
        self.cfIa = 0.  # ia.. incorrectly annotated (there is an overlap with correct gene)
        self.cfTa = 0.  # ta.. true annotation

        # incorrect family annotation
        self.ifCa = 0.  # ca.. correctly annotated
        self.ifIa = 0.  # ia.. incorrectly annotated (there is an overlap with correct gene)
        self.ifTa = 0.  # ta.. true annotation


class ResultHolder():
    def __init__(self):
        self.famInfoDict = {}  # map: gene-family -> FamInfo
        self.taTotal = 0.  # true gene annotation total (bp)

        self.trueGeneNameSet = set()
        self.cfGeneNameSet = set()
        self.ifGeneNameSet = set()

        self.ngAnnotatedTotalBp = 0.
        self.annotatedLengthNoFamList = []  # list of annotation lengths of no family found
        self.cfAnnotatedLenList = []
        self.ifAnnotatedLenList = []

        self.cfOverlapLenList = []  # correct family - list of overlap lengths
        self.ifOverlapLenList = []  # incorrect family - list of overlap lengths
        self.overlapLenTrueList = []  # list of true gene lengths (on reads)


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


def getHmmAnnotationAccuracy(domFile, pfamToGeneFile, scoreThreshold=40, accuracyThreshold=0.8):
    """
        Evaluate the Pfam annotations of the reads.
    """
    # get the names of the other files based on the name of the dom file
    prefix, fileName = os.path.basename(domFile).split('_', 1)
    dirName = os.path.dirname(domFile)
    fqFile = os.path.join(dirName, '%s_%s' % (prefix, 'join.fq.gz'))
    fqProtFile = os.path.join(dirName, '%s_%s' % (prefix, 'join_prot.fna.gz'))
    samFile = os.path.join(dirName, '%s_%s' % (prefix, 'join_gmap.sam.gz'))

    # map: read-name-prot -> seq-prot
    protNameToSeq = fas.fastaFileToDictWholeNames(fqProtFile)

    # map: read-name -> list of hits (lists sorted according to the scores)
    nameToHitList = getReadNameToHitList(domFile)

    # map: read-name -> set of gene-names-true
    nameToTrueGeneSet = {}

    # map: read-name -> list of (rPosOnRead, rPosOnGene, overlapLen, gStrain, gName) (sorted ~ overlap lengths)
    nameToTrueGeneList = {}

    # for each line in the sam-file, store the true read to gene mapping
    for line in gzip.open(samFile):
        line = line.strip()
        if line.startswith('#'):
            continue
        tokens = line.split('\t')

        # the read contains a gene
        if len(tokens) > 11:
            entries = tokens[11].split('|')  # list of gene entries
            readName = tokens[0]
            eList = []
            # get entries as lists
            for e in entries:
                eList.append(e.split(','))
            assert len(eList) >= 1

            # a read maps onto several genes, sort according to the overlap lengths
            if len(eList) > 1:
                eList.sort(key=lambda x: int(x[2]), reverse=True)

            # collect the set of genes
            geneSet = set()
            for e in eList:
                geneSet.add(e[4])
            nameToTrueGeneSet[readName] = geneSet
            nameToTrueGeneList[readName] = eList

    # map: pfam-dom-name -> list of (gene-name, support)
    pfamToGeneList = {}
    for line in open(pfamToGeneFile):
        line = line.strip()
        if line.startswith('#'):
            continue
        tokens = line.split('\t')
        assert len(tokens) >= 2
        if len(tokens) != 2:
            print('Two tokens expected, got: %s' % len(tokens))
        pfamDom = tokens[0]
        geneScore = tokens[1].split(',')

        # list of (gene-name, score)
        geneScoreList = []
        for i in range(len(geneScore) / 2):
            geneScoreList.append((geneScore[2*i], float(geneScore[2*i+1])))
        pfamToGeneList[pfamDom] = geneScoreList

    # compute how good is the HMM annotation (how good is the mapping)
    rh = ResultHolder()

    # go over all reads
    for readName, dna, p, qs in fq.ReadFqGen(fqFile):
        readName = readName[1:]  # strip starting '@'

        # is it HMM annotated?
        topHit = None
        famInfo = None
        if readName in nameToHitList:
            topHit = nameToHitList[readName][0]  # take the hit with the highest score

            # is the hit significant, filter according to the score and accuracy
            if float(topHit[13]) < scoreThreshold or float(topHit[21]) < accuracyThreshold:
                continue
            else:
                # get the fam-info
                famName = topHit[3]
                if famName not in rh.famInfoDict:
                    rh.famInfoDict[famName] = FamInfoHolder(famName)
                famInfo = rh.famInfoDict[famName]

        # does it have any true annotation? (any gene on it?)
        if readName in nameToTrueGeneList:

            # store true annotation total (sum up overlap lengths of all genes on the read)
            for e in nameToTrueGeneList[readName]:
                rh.trueGeneNameSet.add(e[4])  # (rPosOnRead, rPosOnGene, "overlapLen", gStrain, gName)
                overlapLenTrue = int(e[2])
                rh.taTotal += overlapLenTrue
                rh.overlapLenTrueList.append(overlapLenTrue)

            # the hit is on a read containing genes
            if topHit is not None:

                # top hit coordinates within the read
                startOnReadH, overlapLenH, strainH, scoreH, accH = dnaHitInfo(topHit, dna, protNameToSeq[topHit[0]])

                # get the set of genes corresponding to the pfam
                pfamGeneSet = None
                pfamGeneList = pfamToGeneList.get(topHit[3])
                if pfamGeneList is not None:
                    pfamGeneSet = set(map(lambda x: x[0], pfamGeneList))

                # compute overlap correctness between genes and HMM annotation, go over all genes on the read
                for rPosOnReadT, rPosOnGeneT, overlapLenT, gStrainT, gNameT in nameToTrueGeneList[readName]:
                    rPosOnReadT = int(rPosOnReadT)
                    overlapLenT = int(overlapLenT)
                    # gStrainT = int(gStrainT)  # TODO: check also the strain?

                    # get overlap positions
                    overlapStart = max(rPosOnReadT, startOnReadH)
                    overlapEnd = min(rPosOnReadT + overlapLenT - 1, startOnReadH + overlapLenH - 1)

                    # there is an overlap between a gene and an HMM annotation
                    if overlapStart <= overlapEnd:
                        overlap = overlapEnd - overlapStart + 1
                    else:
                        overlap = 0

                    # compute overlap complement length (not matching positions)
                    complement = max(0, rPosOnReadT - startOnReadH) \
                        + max(0, (startOnReadH + overlapLenH - 1) - (rPosOnReadT + overlapLenT - 1))

                    # is the gene and annotated pfam compatible
                    if overlap > 0:
                        if (pfamGeneSet is not None) and (gNameT in pfamGeneSet):
                            famInfo.cfCa += overlap
                            famInfo.cfTa += overlapLenT
                            famInfo.cfIa += complement
                            rh.cfOverlapLenList.append(overlap)
                            rh.cfGeneNameSet.add(gNameT)
                            rh.cfAnnotatedLenList.append(overlapLenH)
                        else:
                            # incompatible gene overlap
                            famInfo.ifCa += overlap
                            famInfo.ifTa += overlapLenT
                            famInfo.ifIa += complement
                            rh.ifOverlapLenList.append(overlap)
                            rh.ifGeneNameSet.add(gNameT)
                            rh.ifAnnotatedLenList.append(overlapLenH)

        # is it HMM annotated, even though there is no gene on it?
        elif topHit is not None:
            # store ia (sum up overlap length of the best scoring hit that is wrong)
            annotatedLengthNoFam = (int(topHit[20]) - int(topHit[19]) + 1) * 3  # (envTo - envFrom + 1) * 3
            famInfo.ngAnnotatedTotalBp += annotatedLengthNoFam
            rh.annotatedLengthNoFamList.append(annotatedLengthNoFam)

    return rh


def mergeResultsHmmAnnotationAccuracy(resultHolderList, outFile=None, scoreThreshold=None, accuracyThreshold=None):
    """

        @param resultHolderList:
        @type resultHolderList: list[ResultHolder]
    """
    # merge the results from all the result holders

    rh = ResultHolder()  # containing all the results merged

    # for all result holders
    for r in resultHolderList:

        # for all PROT families
        for famName, famInfo in r.famInfoDict.iteritems():

            if famName in rh.famInfoDict:
                assert rh.famInfoDict[famName].geneFamily == famInfo.geneFamily
                rh.famInfoDict[famName].ngAnnotatedTotalBp += famInfo.ngAnnotatedTotalBp
                rh.famInfoDict[famName].cfCa += famInfo.cfCa
                rh.famInfoDict[famName].cfIa += famInfo.cfIa
                rh.famInfoDict[famName].cfTa += famInfo.cfTa
                rh.famInfoDict[famName].ifCa += famInfo.ifCa
                rh.famInfoDict[famName].ifIa += famInfo.ifIa
                rh.famInfoDict[famName].ifTa += famInfo.ifTa
            else:
                rh.famInfoDict[famName] = famInfo

        rh.taTotal += r.taTotal

        rh.trueGeneNameSet.update(r.trueGeneNameSet)
        rh.cfGeneNameSet.update(r.cfGeneNameSet)
        rh.ifGeneNameSet.update(r.ifGeneNameSet)

        rh.ngAnnotatedTotalBp += r.ngAnnotatedTotalBp

        rh.annotatedLengthNoFamList.extend(r.annotatedLengthNoFamList)
        rh.cfAnnotatedLenList.extend(r.cfAnnotatedLenList)
        rh.ifAnnotatedLenList.extend(r.ifAnnotatedLenList)

        rh.cfOverlapLenList.extend(r.cfOverlapLenList)
        rh.ifOverlapLenList.extend(r.ifOverlapLenList)
        rh.overlapLenTrueList.extend(r.overlapLenTrueList)

    # get overall statistics

    # correct family annotation
    cfCa = 0.  # ca.. correctly annotated
    cfIa = 0.  # ia.. incorrectly annotated
    cfTa = 0.  # ta.. true annotation

    # incorrect family annotation
    ifCa = 0.  # ca.. correctly annotated
    ifIa = 0.  # ia.. incorrectly annotated
    ifTa = 0.  # ta.. true annotation
    ngAnnotatedTotalBp = 0.  # annotated, but there is nothing there

    # go over all families
    for info in rh.famInfoDict.values():
        cfCa += info.cfCa
        cfIa += info.cfIa
        cfTa += info.cfTa
        ifCa += info.ifCa
        ifIa += info.ifIa
        ifTa += info.ifTa
        ngAnnotatedTotalBp += info.ngAnnotatedTotalBp

    buff = ''
    # print out overall statistics
    if scoreThreshold is not None and accuracyThreshold is not None:
        buff += str('# Score: %s Acc: %s\n' % (scoreThreshold, accuracyThreshold))

    # precision, recall, average overlap length
    buff += str('# CF: P: %s R: %s\n' % (round((cfCa / (cfCa + cfIa)) * 100., 2), round((cfCa / cfTa) * 100., 2)))
    buff += str('# IF: P: %s R: %s\n' % (round((ifCa / (ifCa + ifIa)) * 100., 2), round((ifCa / ifTa) * 100., 2)))

    buff += str('# ---')
    # Portion from true, portion from alignments
    buff += str('# CF: Ttotal: %s Atotal: %s\n' % (round((cfTa / rh.taTotal) * 100., 2), round((cfCa + cfIa) / (cfCa + cfIa + ifCa + ifIa + ngAnnotatedTotalBp) * 100., 2)))
    buff += str('# IF: Ttotal: %s Atotal: %s\n' % (round((ifTa / rh.taTotal) * 100., 2), round((ifCa + ifIa) / (cfCa + cfIa + ifCa + ifIa + ngAnnotatedTotalBp) * 100., 2)))
    buff += str('# NT: Ttotal: %s Atotal: %s\n' % (round(((rh.taTotal - cfTa - ifTa) / rh.taTotal) * 100., 2), round(ngAnnotatedTotalBp / (cfCa + cfIa + ifCa + ifIa + ngAnnotatedTotalBp) * 100., 2)))

    buff += str('# ---\n')
    # overlap lengths
    buff += str('# OverlapLen: CF: %s IF: %s\n' % (round(np.average(rh.cfOverlapLenList), 2), round(np.average(rh.ifOverlapLenList), 2)))

    # average annotated length on a read containing no gene!
    buff += str('# AlignmenLen: T: %s N: %s CF: %s IF: %s\n' % (round(np.average(rh.overlapLenTrueList)),
                                        round(np.average(rh.annotatedLengthNoFamList)),
                                        round(np.average(rh.cfAnnotatedLenList)),
                                        round(np.average(rh.ifAnnotatedLenList))))
    # Number of families
    buff += str('# TotalFamilies: %s\n' % len(rh.famInfoDict))
    buff += str('# Genes: T: %s cf: %s if: %s\n' % (len(rh.trueGeneNameSet), len(rh.cfGeneNameSet), len(rh.ifGeneNameSet)))
    buff += str('# ------------------------------------------------------------------------\n')

    if outFile is None:
        print(buff)
    else:
        out = csv.OutFileBuffer(outFile)
        out.writeText(buff)
        out.writeText('# Detailed results \n')
        out.writeText('# FamName, CF_P, CF_R, IF_P, IF_R, NA, NA_count\n')
        tupleList = []
        for famName, f in rh.famInfoDict.iteritems():
            if com.isclose(0, f.cfCa + f.cfIa):
                a1 = 'NA'
            else:
                a1 = round((f.cfCa / (f.cfCa + f.cfIa)) * 100., 2)

            if com.isclose(0, f.cfTa):
                a2 = 'NA'
            else:
                a2 = round((f.cfCa / f.cfTa) * 100., 2)

            if com.isclose(0, f.ifCa + f.ifIa):
                a3 = 'NA'
            else:
                a3 = round((f.ifCa / (f.ifCa + f.ifIa)) * 100., 2)

            if com.isclose(0, f.ifTa):
                a4 = 'NA'
            else:
                a4 = round((f.ifCa / f.ifTa) * 100., 2)

            if com.isclose(0, f.ngAnnotatedTotalBp):
                a5 = 0.
            else:
                a5 = round((f.ngAnnotatedTotalBp / (f.ngAnnotatedTotalBp + f.cfCa + f.cfIa + f.ifCa + f.ifIa)) * 100., 2)

            tupleList.append((famName, a1, a2, a3, a4, a5, f.ngAnnotatedTotalBp))

        tupleList.sort(key=lambda x: x[5])

        for t in tupleList:
            out.writeText(', '.join(map(lambda x: str(x), t)) + '\n')


def dnaHitInfo(hit, dnaSeq, protSeq=None, translTable=11):
    """
        Given the hit entry from the domtblout file for a prot sequence, return the hit coordinates on the dna sequence.

        @type hit: list[str]
        @type dnaSeq: str

        @param protSeq: if equal to a sequence, verification is performed (no verification if None)

        @return: tuple (startOnRead, overlapLen, strain, score, acc)
    """
    # read hit info
    score = float(hit[13])
    fromEnv = int(hit[19]) - 1
    lenEnv = int(hit[20]) - fromEnv
    # toEnv = int(hit[20]) - 1
    acc = float(hit[21])
    frame = int(hit[0][-1:])
    translSeq = None

    # length of the overlap on the dna sequence
    overlapLen = lenEnv * 3

    # frame shift on the positive strain
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


def partitionReads(sampleDir, scoreThreshold, accuracyThreshold, shuffleRandSeed, pfamPartitionedDir, joinedReads=True,
                   considerSam=True):
    """
        Partitioning reads into the individual Pfam-domains.
    """
    try:
        strainDirList = map(lambda x: os.path.join(sampleDir, x), os.listdir(sampleDir))

        samplePartDir = os.path.join(sampleDir, pfamPartitionedDir)
        if not os.path.isdir(samplePartDir):
            os.mkdir(samplePartDir)

        # for each pfam-dom
        fqOutDict = {}
        fqProtOutDict = {}
        fqSamOutDict = {}
        fqDomOutDict = {}

        for strainDir in strainDirList:
            strainAcc = os.path.basename(strainDir)
            if strainAcc == 'sample_partitioned':  # skip the directory containing the partitioned data
                continue
            if os.path.isdir(strainDir):
                # find the dom file! do it for each dom file you find
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
                                        if len(line.split('\t')) == 11:  # lines with only 11 entries will be padded with * to 12
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

                            # go over all reads !!!
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
                                    startOnRead, overlapLen, strain = dnaHitInfo(topHit, dna, protSeq)[:3]

                                    fqDomOutDict[famName].append('\t'.join(topHit) + '\t%s\t%s\t%s' % (startOnRead, overlapLen, strain))

        if joinedReads:
            ident = 'join'
        else:
            ident = 'pair'

        # for each pfam-dom, store into a file
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

            # write fq
            fqOut = fq.WriteFq(fqOutO)
            for e in fqContentList:
                fqOut.writeFqEntry('@' + e[0], e[1], e[2], e[3])
            fqOut.close()

            # write prot
            fqProtOut = fq.WriteFq(fqProtOutO)
            fqProtOut.write('\n'.join(map(lambda x: '>%s\n%s' % (x[0], x[1]), fqProtOutDict[famName])) + '\n')
            fqProtOut.close()

            # write sam
            if considerSam:
                fqSamOut = fq.WriteFq(fqSamOutO)
                fqSamOut.write('\n'.join(fqSamOutDict[famName]) + '\n')
                fqSamOut.close()

            # write dom
            fqDomOut = fq.WriteFq(fqDomOutO)
            fqDomOut.write('\n'.join(fqDomOutDict[famName]) + '\n')
            fqDomOut.close()

            # shuffle file entries
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


def getStatPartitionedReads(sampleDir, pfamPartitionedDir, pfamPartitionedStatFile):

    samplePartDir = os.path.join(sampleDir, pfamPartitionedDir)
    lineList = []

    for f in os.listdir(samplePartDir):
        if f.startswith('r_') and f.endswith('join_prot.domtblout.gz'):
            domFile = os.path.join(samplePartDir, f)
            fTokens = f.split('_')
            gsamFile = os.path.join(samplePartDir, '_'.join(fTokens[:-2]) + '_join_gmap.sam.gz')
            assert os.path.isfile(domFile) and os.path.isfile(gsamFile)
            domEntryCount = 0
            domName = set()
            domAvgScore = 0.
            for line in gzip.open(domFile):
                line = line.strip()
                if not line.startswith('#'):
                    tokens = line.split('\t', 14)
                    if len(tokens) > 14:
                        domEntryCount += 1
                        domName.add(tokens[3])  # famName
                        domAvgScore += float(tokens[13])  # score
            domAvgScore /= float(domEntryCount)

            samEntryCount = 0
            strainSet = set()
            geneSet = set()
            for line in gzip.open(gsamFile):
                line = line.strip()
                if not line.startswith('#'):
                    tokens = line.rsplit('\t', 2)
                    if len(tokens) > 2:
                        samEntryCount += 1
                        strainSet.add(tokens[-1])
                        geneSet.add(tokens[-2].rsplit(',', 1)[-1])

            assert len(domName) == 1  # there should be only one domain
            domName = domName.pop()
            assert domEntryCount == samEntryCount
            assert len(strainSet) > 0
            assert len(geneSet) > 0
            geneList = list(geneSet)  # '*' means no gene
            geneList.sort()
            lineList.append((domName, comh.getGeneNameToFileName(domName)[:-4], domEntryCount, len(strainSet), len(geneSet), round(domAvgScore, 2), '|'.join(geneList)))

    lineList.sort(key=lambda x: x[2], reverse=True)

    out = csv.OutFileBuffer(os.path.join(sampleDir, pfamPartitionedStatFile))
    out.writeText('# famName, famFileEnc, readCount, strainCount, genesCount, avgScore, geneList\n')

    for entry in lineList:
        out.writeText('\t'.join(map(lambda x: str(x), entry)) + '\n')
    out.close()


def testPartition():
    partitionReads('/home/igregor/Documents/work/hsim/562/samples/0',
                   comh.SAMPLES_PFAM_EVAN_MIN_SCORE,
                   comh.SAMPLES_PFAM_EVAN_MIN_ACCURACY,
                   comh.SAMPLES_SHUFFLE_RAND_SEED,
                   comh.SAMPLES_PFAM_PARTITIONED_DIR)


def testPartitionStat():
    getStatPartitionedReads('/home/igregor/Documents/work/hsim/562/samples/0',
                            comh.SAMPLES_PFAM_PARTITIONED_DIR,
                            comh.SAMPLES_PFAM_PARTITIONED_STAT_FILE)


    # avg score? if a gene was located, how good the match was? (score / length)?
    # perfect scores would mean that the reference is too close !!!

    # design measures (precision, recall analoguous)!
    # recall = there is a gene, how good was it recalled/recovered

    # base the reading frame on the best score, look whether there is another match to the using the same
    # reading frame or the reverse complement?


def testAnnotationAccuracy():
    # _main()
    # original sample NZ_AEZU00000000 !!!  next: NZ_AIGN00000000
    domFile = '/home/igregor/Documents/work/hsim/562/samples/1/NZ_AEZU00000000/0_join_prot.domtblout.gz'
    pfamToGeneFile = '/home/igregor/Documents/work/hsim/562/fasta_pull_genes_pfam_to_gene_name.csv'
    accuracyThreshold = 0.8
    accuracyList = [0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99]
    scoreList = [0, 30, 40, 50, 60, 70, 80]
    # scoreList = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]

    # for scoreThreshold in scoreList:
    #     for accuracyThreshold in accuracyList:
    #         getHmmAnnotationAccuracy(domFile, pfamToGeneFile, scoreThreshold, accuracyThreshold)
    #
    mergeResultsHmmAnnotationAccuracy([getHmmAnnotationAccuracy(domFile, pfamToGeneFile, 80, 0.99),
                                       getHmmAnnotationAccuracy(domFile, pfamToGeneFile, 80, 0.99)])

    # name = 'name_5'
    # score = 123.0
    # fromEnv = 3
    # toEnv = 3
    # dnaSeq = 'CTATCGACAGCTGTCTA'
    # protSeq = 'RQLSI'
    # dnaHitInfo([name,1,2,3,4,5,6,7,8,9,10,11,12,score,14,15,16,17,18,fromEnv,toEnv,21], dnaSeq, protSeq)

    # print 'targetName', targetName
    # print 'accession1', accession1
    # print 'tlen', tlen
    # print 'queryName', queryName
    # print 'accession2', accession2
    # print 'qlen', qlen
    # print 'eValue', eValue
    # print 'score1', score1
    # print 'bias', bias
    # print 'cEvalue', cEvalue
    # print 'iEvalue', iEvalue
    # print 'score2', score2
    # print 'bias', bias
    # print 'from1', from1
    # print 'to1', to1
    # print 'from2', from2
    # print 'to2', to2
    # print 'from3', from3
    # print 'to3', to3
    # print 'acc', acc
    # print 'descriptionOfTarget', descriptionOfTarget
    # print 'rest', r1, r2

if __name__ == "__main__":
    _main()
    # testAnnotationAccuracy()
    # testPartition()
    # testPartitionStat()