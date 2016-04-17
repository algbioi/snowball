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

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.

    Contains constants and common functionality used in this package. Module version 1.2
"""
import os
import sys
import re
import tarfile
import multiprocessing as mp
from algbioi.com import parallel
import platform

# list of species to be processed
SPECIES_LIST = ["562"]  # Escherichia coli, species
# SPECIES_LIST = ["1166683"]  # Rhizobia

# Art simulator setting (only the first val used)
ART_READ_LEN = [150, 150]
ART_INSERT_SIZE = [225, 5000]  # 150, 180, 450  # try 180 and 5000
ART_INSERT_SD = [23, 500]  # "10 %"
ART_MIN_SEQ_LEN = [200, 5000]  # this should be handled by art !!!
ART_QS_MAX = [60, 60]  # maximum quality score
ART_Q_PROFILE = [(None, None), (None, None)]  # the default profiles are taken !

# maximum number of CPUs to be used
MAX_PROC = min(12, mp.cpu_count())

TRANSLATION_TABLE = 11

# Data locations
if sys.platform == 'darwin':
    # REFERENCE_DIR_ROOT = '/Volumes/VerbatimSSD/work/hsim01'  # external storage !!!
    REFERENCE_DIR_ROOT = '/Users/ivan/Documents/nobackup/hsim01'  # local disk !!!
    NCBI_TAXONOMY_FILE = '/Users/ivan/Documents/work/binning/taxonomy/20140916/ncbitax_sqlite.db'
    HMM_PROFILE_DIR = '/Volumes/VerbatimSSD/work/pfam'
else:
    if platform.dist()[0] == 'Ubuntu':
        REFERENCE_DIR_ROOT = '/home/igregor/Documents/work/hsim'
        NCBI_TAXONOMY_FILE = '/home/igregor/Documents/work/taxonomy/20140916/ncbitax_sqlite.db'
        HMM_PROFILE_DIR = '/home/igregor/Documents/work/db/pfamV27'
    else:
        assert sys.platform == 'linux2'
        REFERENCE_DIR_ROOT = '/net/metagenomics/projects/PPSmg/hsim/hsim01'
        NCBI_TAXONOMY_FILE = '/net/metagenomics/projects/PPSmg/taxonomy/20140916/ncbitax_sqlite.db'
        HMM_PROFILE_DIR = '/net/metagenomics/projects/PPSmg/database/pfamV27/nobackup/pub/databases/Pfam/releases/Pfam27.0'

# Directory names
FASTA_GENOMES_DIR_NAME = 'fasta_genomes'
FASTA_GENOMES_DRAFT_DIR_NAME = 'fasta_draft_genomes'

FASTA_GENOMES_GENES_DIR_NAME = 'fasta_genomes_genes'
FASTA_GENOMES_DRAFT_GENES_DIR_NAME = 'fasta_draft_genomes_genes'
FASTA_PULL_GENES_DIR_NAME = 'fasta_pull_genes'

FASTA_PULL_GENES_PHYLO_DIR_NAME = 'fasta_pull_genes_phylo'
FASTA_PULL_GENES_PHYLO_ALIGN_DIR_NAME = 'fasta_pull_genes_phylo_align'
FASTA_PULL_GENES_PHYLO_PROT_DIR_NAME = 'fasta_pull_genes_phylo_prot'
FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME = 'fasta_pull_genes_pfam_to_gene_name.csv'
FASTA_PULL_GENES_PHYLO_PFAM_TO_GENE_NAME_AMBIGUOUS = 'fasta_pull_genes_pfam_to_gene_name_ambiguous.csv'
FASTA_PULL_GENES_PHYLO_GENE_NAME_TO_PFAM = 'fasta_pull_genes_gene_name_to_pfam.csv'


CLUSTER_PHYLO_DIR = 'cluster_phylo'
CLUSTER_PHYLO_STRAINS_FNA = 'genes_concat_align.fna'
CLUSTER_MOTHUR_METHOD = 'furthest'  # only this one defined
CLUSTER_MOTHUR_CWD = 'mothur_cwd'

SAMPLES_DIR = 'samples'
SAMPLES_DEF_FILE = 'samples_definitions.txt'  # individual sample definitions
SAMPLES_ERROR_PROFILE = 'samples_error_profile.csv'
SAMPLES_ERROR_QS_CUTOFF = 'samples_error_qs_cutoffs.csv'
SAMPLES_JOIN_ERROR_PROFILE = 'samples_join_error_profile.csv'
SAMPLES_JOIN_PFAM_MAP_QUALITY = 'samples_join_pfam_map_quality.csv'
SAMPLES_PFAM_PARTITIONED_DIR = 'sample_partitioned'
SAMPLES_PFAM_PARTITIONED_STAT_FILE = 'sample_partitioned_stat.csv'
SAMPLES_SHUFFLE_RAND_SEED = 1

SAMPLES_DEF_STRAIN_COUNT_LIST = (3,)  # (3, 4, 5, 6, 7, 8, 9)
SAMPLES_DEF_IDENTITY_CUTOFF = 0.002
SAMPLES_DEF_MAX_SAMPLES_PER_CLUSTER = 1  # 3
SAMPLES_DEF_RAND_SEED = 1
SAMPLES_DEF_LOGNORM_MEAN = 1
SAMPLES_DEF_LOGNORM_SD = 2
SAMPLES_DEF_MIN_COVERAGE = 1
SAMPLES_DEF_MAX_COVERAGE = 50

# for read filtering .. not used?
SAMPLES_READ_TRIM_CUTOFFS = [0.25, 0.25]  # for each library: levels at which the QS cutoffs will be considered
SAMPLES_READ_TRIM_REMAIN = [0.5, 0.5]  # this continuous part of a read must remain, else it's thrown away

# joining pair end
SAMPLES_PAIRED_END_JOIN_MIN_OVERLAP = 0.05
SAMPLES_PAIRED_END_JOIN_MIN_OVERLAP_IDENTITY = 0.9

SAMPLES_PFAM_EVAN_MIN_SCORE = 40
SAMPLES_PFAM_EVAN_MIN_ACCURACY = 0.6
HMM_PROFILE_FILE = 'Pfam-A_and_Amphora2.hmm'
HMM_PROFILE_FILES_PARTITIONED_DIR = 'pfam_a_and_amphora2'

# Assembly parameters
ASSEMBLY_CONSIDER_PROT_COMP = (False,)
ASSEMBLY_ONLY_POVERLAP = (False,)
ASSEMBLY_POVERLAP = (0.8,)  # 0.7 - 0.8
ASSEMBLY_OVERLAP_LEN = (0.5,)  # 0.25 - 0.8
ASSEMBLY_OVERLAP_ANNOT_LEN = (0.2,)
ASSEMBLY_STOP_OVERLAP_MISMATCH = (0.1,)
ASSEMBLY_MAX_LOOPS = 1
ASSEMBLY_SUPER_READ_EVAL_INIT = 'super_read_init_eval.txt'
ASSEMBLY_SUPER_READ_STAT_INIT = 'super_read_init_stat.txt'
ASSEMBLY_SUPER_READ_EVAL_INIT_SAT = 'sat_super_read_init_eval.txt'
ASSEMBLY_SUPER_READ_EVAL_CLEN = 'super_read_clen.txt'
ASSEMBLY_SUPER_READ_EVAL_CLEN_SAT = 'sat_super_read_clen.txt'
ASSEMBLY_SUPER_READ_EVAL_RCOV = 'super_read_rcov.txt'
ASSEMBLY_SUPER_READ_EVAL_RCOV_SAT = 'sat_super_read_rcov.txt'


# Binary locations
if sys.platform == 'darwin':
    MUSCLE_BINARY = '/Users/ivan/Documents/work/tools/muscle/muscle3.8.31_i86darwin64'
    MOTHUR_BINARY = '/Users/ivan/Documents/work/tools/mothur/mothur/mothur'
    ART_ILLUMINA_BINARY = '/Users/ivan/Documents/work/tools/art/art_bin_ChocolateCherriesOSX/art_illumina'
    HMMER_BINARY = '/Users/ivan/Documents/work/tools/hmmer/hmmer-3.0-macosx-intel/binaries'
    SAT_INSTALL_DIR = '/Users/ivan/Documents/work/tools/sat/SAT-Assembler'
    BWA_INSTALL_DIR = '/Users/ivan/Documents/work/tools/Bowtie/bowtie2-2.2.5'
else:
    if platform.dist()[0] == 'Ubuntu':
        MUSCLE_BINARY = '/home/igregor/Documents/work/tools/muscle/muscle3.8.31_i86linux64'
        MOTHUR_BINARY = '/home/igregor/Documents/work/tools/mothur/mothur'
        ART_ILLUMINA_BINARY = '/home/igregor/Documents/work/tools/art/art_bin_ChocolateCherriesLinux/art_illumina'
        HMMER_BINARY = '/home/igregor/Documents/work/tools/hmmer-3.0/binaries'
        SAT_INSTALL_DIR = '/home/igregor/Documents/work/tools/sat/SAT-Assembler'
        BWA_INSTALL_DIR = '/home/igregor/Documents/work/tools/bowtie/bowtie2-2.2.5'
    else:
        assert sys.platform == 'linux2'
        MUSCLE_BINARY = '/net/metagenomics/projects/PPSmg/hsim/muscle3.8.31_i86linux64'
        MOTHUR_BINARY = '/net/metagenomics/projects/PPSmg/tools/mothur/mothur_1_333/mothur'
        ART_ILLUMINA_BINARY = '/net/metagenomics/projects/PPSmg/tools/art/art_bin_ChocolateCherriesLinux/art_illumina'
        HMMER_BINARY = '/net/metagenomics/projects/PPSmg/tools/hmmer-3.0/binaries'
        SAT_INSTALL_DIR = '/net/metagenomics/projects/PPSmg/tools/sat/SAT-Assembler'
        BWA_INSTALL_DIR = '/net/metagenomics/projects/PPSmg/tools/Bowtie/bowtie2-2.2.5'

# Common functionality

def getGeneNameToFileName(geneName):
    """
        From a gene name get a FASTA file name.
        All letters lowered, chars [^a-z0-9] are replaced by underscore)
    """
    fileName = re.sub(r"[^a-z0-9]", "_", geneName.lower())
    changesNum = 0
    for i, j in zip(list(geneName), list(fileName)):
        if i != j:
            changesNum += 1
    fileName = fileName + '_' + str(changesNum) + '.fna'
    return fileName


def getSeqNameToDict(seqName):
    """
        Converts a sequence name (e.g. >key1:val1;key2:val2) to a dictionary.

        type seqName: str
        rtype: dict
    """
    d = {}
    for e in seqName.lstrip('>').split(';'):
        k, v = e.split(':')
        assert k not in d
        d[k] = v
    return d


def getAlignments(inDir, outDir, reportFailedCmd=True, listOfInFileNames=None):
    """
        Build a multiple sequence alignments for FASTA files in the input directory.

        @param inDir: directory containing input fasta files
        @param outDir: directory for the output aligned fasta files
        @param listOfInFileNames: list of file names from the input directory for which the alignments will be computed
        (default ~ None means for all files in the directory)
        @return: a list of failed commands or None
    """
    assert os.path.isfile(MUSCLE_BINARY), 'Binary file does not exist: %s' % MUSCLE_BINARY
    # compute alignments for all files if the list of files is not defined
    if listOfInFileNames is None:
        listOfInFileNames = os.listdir(inDir)
    # define a task list
    taskList = []
    assert os.path.isdir(outDir)
    for fileName in listOfInFileNames:
        assert os.path.isfile(os.path.join(inDir, fileName))
        cmd = '%s -in %s -out %s' % (MUSCLE_BINARY, os.path.join(inDir, fileName), os.path.join(outDir, fileName))
        taskList.append(parallel.TaskCmd(cmd, outDir))
    # run tasks in parallel
    failedCmd = parallel.runCmdParallel(taskList, MAX_PROC)
    if reportFailedCmd:
        parallel.reportFailedCmd(failedCmd)
    return failedCmd


def getPhyloGenesFilePath(specDir):
    """
        Given the species directory, returns a file containing gene names for the phylogeny reconstruction / clustering.
    """
    phyloFile = None
    for f in os.listdir(specDir):
        if f.startswith('phylo_genes_'):
            assert phyloFile is None
            phyloFile = os.path.join(specDir, f)
    assert os.path.isfile(phyloFile)
    return phyloFile


def isSpeciesDirectory(dirName, taxonomy):
    """
        @type taxonomy: TaxonomyNcbi
    """
    try:
        if taxonomy.getRank(int(os.path.basename(dirName))) == 'species':
            return True
    except ValueError:
        pass
    return False


def extract(compressedFile):
    """
        Decompress a file with a ".tgz" ending.
        A file won't be decompressed in the case it has already been decompressed.
    """
    assert compressedFile.endswith('.tgz'), 'Not supported extension: %s' % compressedFile
    dstDir = os.path.normpath(compressedFile.replace('.tgz', ''))
    if not os.path.isdir(dstDir):
        tar = tarfile.open(compressedFile, 'r:gz')
        for item in tar:
            tar.extract(item, dstDir)


class SampleDef(object):
    def __init__(self, defFile):
        """
            @param defFile: SAMPLES_DEF_FILE
        """
        self.idToStrainList = {}
        for line in open(defFile):
            line = line.strip()
            if len(line) > 0:
                tokens = line.split('\t', 2)
                sId = int(tokens[0])
                strains = tokens[1].split(',')
                self.idToStrainList[sId] = strains

    def sIdToStrainListLen(self, sId):
        """
            @return: of how many strains a sample consist
        """
        strains = self.idToStrainList.get(int(sId))
        if strains is not None:
            return len(strains)
        else:
            return None
