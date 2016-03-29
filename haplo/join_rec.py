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

    Implements joining of the reads for the "Snowball algorithm".
"""

# TODO: implement joining of the reads

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


def getCodingTable(translTable=11):
    nucL = ['A', 'T', 'G', 'C']
    buff = ''
    aToN = {}
    for p1 in nucL:
        for p2 in nucL:
            for p3 in nucL:
                dna = p1 + p2 + p3
                prot = str(Seq(dna, generic_dna).translate(table=translTable))

                # try:
                #     prot = str(Seq(dna + 'TAG', generic_dna).translate(table=translTable, stop_symbol='', cds=True))
                #     print 'Start:' + dna + ' : ' +  prot
                # except Exception as e:
                #     print e.message
                    # pass
                if prot not in aToN:
                    aToN[prot] = []
                aToN[prot].append(dna)

                # buff += '%s : %s\n' % (dna, prot)
    for k, v in aToN.iteritems():
        buff += k + ': ' + ','.join(v) + '\n'

    return buff

