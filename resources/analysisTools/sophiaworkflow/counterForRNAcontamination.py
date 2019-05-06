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
from sys import exit as sysexit

inputPath = argv[1]
eventDict = dict()
with open(inputPath) as inputHandle:
    for line in inputHandle:
        lineChunks = line.rstrip().split('\t')
        direction = int(lineChunks[3])
        lineIndex = int(lineChunks[4])
        gene = lineChunks[-1].split('_')[0]
        if lineIndex not in eventDict:
            eventDict[lineIndex] = [set(), set()]
        eventDict[lineIndex][direction].add(gene)

resultDict = dict()
for key in eventDict:
    event = eventDict[key]
    if len(event[0]) > 0 and len(event[1]) > 0:
        intersectionGenes = event[0].intersection(event[1])
        if len(intersectionGenes) > 0:
            for gene in intersectionGenes:
                if gene not in resultDict:
                    resultDict[gene] = 0
                resultDict[gene] += 1
for gene in resultDict:
    print(gene, resultDict[gene], sep='\t')
print("#done")
sysexit(0)
