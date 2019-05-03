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

broadProd = list(product(range(51), range(51)))
narrowProd = list(product(range(10), range(10)))

lineIndex = 0
with open(preFilteredBedpe) as inputHandle:
    for line in inputHandle:
        if line[0] != '#':
            lineChunks = line.rstrip().split('\t')
            eventType = lineChunks[8]
            if lineChunks[0] == lineChunks[3] and eventType in {"TRA", "DEL"}:
                minPos = min(int(lineChunks[1]), int(lineChunks[4]))
                maxPos = max(int(lineChunks[1]), int(lineChunks[4]))
                eventScore = int(lineChunks[9])
                for iterPair in broadProd:
                    if (minPos - iterPair[0]) < maxPos - 1 + iterPair[1]:
                        print(lineChunks[0], minPos - iterPair[0], maxPos - 1 + iterPair[1], lineIndex, sep='\t')
            lineIndex += 1
