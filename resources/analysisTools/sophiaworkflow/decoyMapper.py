#!/usr/bin/env python
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

from sys import argv
from bisect import bisect_left

decoyRangeRef = argv[1]
inputFile = argv[2]

decoyStartPositions = []
decoyEndPositions = []
decoyData = []


def sourceQualCheck(sa):
    if len(sa) == 1:
        return False
    if "-" in sa or "?" in sa:
        return False
    support, secondarySupport, _ = sa.split('(')[1].split(',')
    support = int(support)
    if support < 2:
        return False
    secondarySupport = int(secondarySupport)
    if secondarySupport < 2:
        return False
    return True


with open(decoyRangeRef) as f:
    for line in f:
        decoyStartPosition, decoyEndPosition, mappedChr, mappedStartPos, mappedEndPos = line.rstrip().split('\t')
        decoyStartPositions.append(int(decoyStartPosition))
        decoyEndPositions.append(int(decoyEndPosition))
        decoyData.append([mappedChr, mappedStartPos, mappedEndPos])


def decoyCheck(inputStr):
    if "hs" in inputStr:
        return True
    if "GL" in inputStr:
        return True
    if "Y" in inputStr:
        return True
    if "NC" in inputStr:
        return True
    return False


with open(inputFile) as f:
    for line in f:
        lineChunks = line.rstrip().split('\t')
        oldChr1 = lineChunks[0]
        oldStart1 = lineChunks[1]
        oldEnd1 = lineChunks[2]
        oldChr2 = lineChunks[3]
        oldStart2 = lineChunks[4]
        oldEnd2 = lineChunks[5]
        source1 = lineChunks[16]
        source2 = lineChunks[17]
        if "hs37d5" in {lineChunks[0], lineChunks[3]}:
            if lineChunks[0] == "hs37d5":
                index = bisect_left(decoyStartPositions, int(lineChunks[1]))
                if index > 0:
                    newChr, newPos, newPos2 = decoyData[index - 1]
                    if newChr != "hs37d5":
                        if lineChunks[3] == newChr:
                            newSize = min(abs(int(lineChunks[1]) - int(newPos)), abs(int(lineChunks[2]) - int(newPos2)))
                            if newSize < 1e4:
                                if "-" in source1 or "-" in source2 or "?" in source1 or "?" in source2:
                                    lineChunks[9] = "1"
                        if sourceQualCheck(source2) and (int(newPos2) - int(newPos)) < 21000:
                            lineChunks[0] = newChr
                            lineChunks[1] = newPos
                            lineChunks[2] = newPos2
                            if lineChunks[3] == newChr:
                                lineChunks[10] = str(newSize)
                            else:
                                lineChunks[8] = "TRA"
                                lineChunks[10] = "NA"
            if lineChunks[3] == "hs37d5":
                index = bisect_left(decoyStartPositions, int(lineChunks[4]))
                if index > 0:
                    newChr, newPos, newPos2 = decoyData[index - 1]
                    if newChr != "hs37d5":
                        if lineChunks[0] == newChr:
                            newSize = min(abs(int(lineChunks[1]) - int(newPos)), abs(int(lineChunks[2]) - int(newPos2)))
                            if newSize < 1e4:
                                if "-" in source1 or "-" in source2 or "?" in source1 or "?" in source2:
                                    lineChunks[9] = "1"
                        if sourceQualCheck(source1) and (int(newPos2) - int(newPos)) < 21000:
                            lineChunks[3] = newChr
                            lineChunks[4] = newPos
                            lineChunks[5] = newPos2
                            if lineChunks[0] == newChr:
                                lineChunks[10] = str(newSize)
                            else:
                                lineChunks[8] = "TRA"
                                lineChunks[10] = "NA"
        lineChunks.extend([oldChr1, oldStart1, oldEnd1, oldChr2, oldStart2, oldEnd2])
        print(*lineChunks, sep='\t')
