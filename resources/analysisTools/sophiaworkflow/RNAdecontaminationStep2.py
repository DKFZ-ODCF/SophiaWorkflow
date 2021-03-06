#!/usr/bin/env python
#
# Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


from sys import argv, stderr
from itertools import product

preFilteredBedpe = argv[1]
lineIndicesFile = argv[2]


def getIntronExonIndex(geneRawName):
    geneChunks = geneRawName.split('_')
    geneName = "."
    if geneChunks[0] != ".":
        geneName = '_'.join(geneChunks[:-1])
    index = 1
    intronStatus = "intron" in geneChunks[-1] or geneName == "."
    if geneChunks[-1] != ".":
        if "intron" in geneChunks[-1]:
            index = int(geneChunks[-1].split('intron')[-1])
        elif "exon" in geneChunks[-1]:
            index = int(geneChunks[-1].split('exon')[-1])
    return [geneName, index, intronStatus]


lineIndices = set()
with open(lineIndicesFile) as f:
    for line in f:
        lineIndex = int(line.rstrip())
        lineIndices.add(lineIndex)

lineIndex = -1
with open(preFilteredBedpe) as inputHandle:
    for line in inputHandle:
        if line[0] != '#':
            lineIndex += 1
            skipLine = False
            if lineIndex not in sorted(list(lineIndices)):
                lineChunks = line.rstrip().split('\t')
                eventType = lineChunks[8]
                eventScore = int(lineChunks[9])
                if lineChunks[0] == lineChunks[3] and lineChunks[11] != "INV":
                    gene1Raw = lineChunks[20]
                    gene2Raw = lineChunks[30]
                    gene1ListPre = [getIntronExonIndex(x) for x in [y.split(';')[0] for y in gene1Raw.split(',')] if '_' in x]
                    gene2ListPre = [getIntronExonIndex(x) for x in [y.split(';')[0] for y in gene2Raw.split(',')] if '_' in x]
                    for eventPair in product(gene1ListPre, gene2ListPre):
                        if eventPair[0][0] == eventPair[1][0]:
                            if eventPair[0][0] != ".":
                                if (eventPair[0][2] or eventPair[1][2]):
                                    if abs(eventPair[0][1] - eventPair[1][1]) < 2:
                                        skipLine = True
                                        break
                                else:
                                    if eventType == "DEL" and abs(eventPair[0][1] - eventPair[1][1]) == 1:
                                        skipLine = True
                                        break
                if not skipLine:
                    print(line.rstrip())
        else:
            print(line.rstrip())
