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

from sys import argv

sophiaOutput = argv[1]
defaultReadLength = int(argv[2])


class SuppAlignment:
    def __init__(self, saStr):
        self.invalid = saStr == "_"
        if not self.invalid:
            self.suspiciousMapping = '?' in saStr
            self.overcorrectedProperPairing = '#' in saStr
            saChunks = saStr.replace('#', '').replace('?', '').rstrip(')').split('(')
            posPart = saChunks[0]
            eviPart = saChunks[1]
            self.rightSplit = posPart[0] == '|'
            self.leftSplit = posPart[-1] == '|'
            self.inverted = '_INV' in posPart
            posPart = posPart.lstrip('|').rstrip('|').rstrip('_INV')
            self.chr, pos = posPart.split(':')
            if '-' in pos:
                self.fuzzy = True
                self.startPos, self.endPos = pos.split('-')
                self.startPos = int(self.startPos)
                self.endPos = int(self.endPos)
            else:
                self.fuzzy = False
                self.startPos = int(pos)
                self.endPos = self.startPos
            self.softClip, self.hardClip, mateEvidence = eviPart.split(',')
            self.softClip = int(self.softClip)
            self.hardClip = int(self.hardClip)
            self.mateSupport, self.expectedMateSupport = mateEvidence.split('/')
            self.mateSupport = int(self.mateSupport)
            self.expectedMateSupport = int(self.expectedMateSupport)

    def absolutelySuperiorTo(self, sa2):
        return not self.invalid and not sa2.invalid and self.softClip > sa2.softClip and self.hardClip > sa2.hardClip and self.mateSupport > sa2.mateSupport

    def isLowQualFuzzy(self):
        res = (not self.invalid and
               self.fuzzy and
               self.softClip < 2 and
               self.hardClip < 2 and
               not (self.softClip > 0 and self.hardClip > 0)
               )
        return res

    def __eq__(self, sa2):
        equality = (not self.invalid and
                    not sa2.invalid and
                    self.fuzzy == sa2.fuzzy and
                    self.leftSplit == sa2.leftSplit and
                    self.rightSplit == sa2.rightSplit and
                    self.startPos == sa2.startPos and
                    self.endPos == sa2.endPos and
                    self.chr == sa2.chr and
                    self.softClip == sa2.softClip and
                    self.hardClip == sa2.hardClip and
                    self.mateSupport == sa2.mateSupport and
                    self.expectedMateSupport == sa2.expectedMateSupport and
                    self.inverted == sa2.inverted)
        return equality

    def fuzzySaOverlap(self, nonFuzzySa, offset):
        if nonFuzzySa.invalid:
            return False
        if self.chr != nonFuzzySa.chr:
            return False
        if self.startPos <= nonFuzzySa.startPos <= self.endPos:
            return True
        if self.startPos - offset <= nonFuzzySa.startPos <= self.startPos:
            return True
        if self.endPos <= nonFuzzySa.startPos <= self.endPos + offset:
            return True
        return False


def absoluteSuperiority(chr_1, chr_2, pos_1, pos_2, sa_1, sa_2, defaultReadLength):
    if chr_1 != chr_2:
        return 0
    if pos_1 == pos_2:
        return 0
    if sa_1.invalid or sa_2.invalid:
        return 0
    if abs(pos_1 - pos_2) > defaultReadLength:
        return 0
    if sa_1.fuzzy or sa_2.fuzzy:
        return 0
    if sa_1.startPos != sa_2.startPos:
        return 0
    if sa_1.leftSplit:
        if not sa_2.leftSplit:
            return 0
        if pos_2 > pos_1:
            if sa_1.absolutelySuperiorTo(sa_2):
                return 1
            else:
                return 0
        elif pos_1 > pos_2:
            if sa_2.absolutelySuperiorTo(sa_1):
                return 2
            else:
                return 0
    elif sa_1.rightSplit:
        if not sa_2.rightSplit:
            return 0
        if pos_2 > pos_1:
            if sa_2.absolutelySuperiorTo(sa_1):
                return 2
            else:
                return 0
        elif pos_1 > pos_2:
            if sa_1.absolutelySuperiorTo(sa_2):
                return 1
            else:
                return 0


def fuzzySaOverlap(fuzzySa, nonFuzzySa, offset):
    if nonFuzzySa == "_":
        return False
    fuzzySaChr = fuzzySa.split(':')[0].lstrip('|')
    nonFuzzySaChr = nonFuzzySa.split(':')[0].lstrip('|')
    if nonFuzzySaChr != fuzzySaChr:
        return False
    nonFuzzyPos = int(nonFuzzySa.split(':')[1].split('(')[0].split('_')[0].rstrip('|'))
    fuzzySaChunks = fuzzySa.split(':')[1].split('-')
    fuzzyStartPos = int(fuzzySaChunks[0])
    fuzzyEndPos = int(fuzzySaChunks[1].split('(')[0].split('_')[0].rstrip('|'))
    if fuzzyStartPos <= nonFuzzyPos <= fuzzyEndPos:
        return True
    if fuzzyStartPos - offset <= nonFuzzyPos <= fuzzyStartPos:
        return True
    if fuzzyEndPos <= nonFuzzyPos <= fuzzyEndPos + offset:
        return True
    return False


def dedupFuzzy(resultsWithDups, dedupOffset):
    markedForRemoval = [False for _ in range(len(resultsWithDups))]
    for i in range(len(resultsWithDups)):
        lineChunksI = resultsWithDups[i]
        sa1i = SuppAlignment(lineChunksI[16])
        sa2i = SuppAlignment(lineChunksI[17])
        for j in range(len(resultsWithDups)):
            if i == j:
                continue
            lineChunksJ = resultsWithDups[j]
            sa1j = SuppAlignment(lineChunksJ[16])
            sa2j = SuppAlignment(lineChunksJ[17])
            if sa1i == sa1j:
                if sa2i.isLowQualFuzzy() and not sa2j.isLowQualFuzzy():
                    if markedForRemoval[i]:
                        continue
                    if sa2i.fuzzySaOverlap(sa2j, dedupOffset):
                        markedForRemoval[i] = True
                elif sa2j.isLowQualFuzzy() and not sa2i.isLowQualFuzzy():
                    if markedForRemoval[j]:
                        continue
                    if sa2j.fuzzySaOverlap(sa2i, dedupOffset):
                        markedForRemoval[j] = True
            elif sa2i == sa2j:
                if sa1i.isLowQualFuzzy() and not sa1j.isLowQualFuzzy():
                    if markedForRemoval[i]:
                        continue
                    if sa1i.fuzzySaOverlap(sa1j, dedupOffset):
                        markedForRemoval[i] = True
                elif sa1j.isLowQualFuzzy() and not sa1i.isLowQualFuzzy():
                    if markedForRemoval[j]:
                        continue
                    if sa1j.fuzzySaOverlap(sa1i, dedupOffset):
                        markedForRemoval[j] = True
            elif sa1i == sa2j:
                if sa2i.isLowQualFuzzy() and not sa1j.isLowQualFuzzy():
                    if markedForRemoval[i]:
                        continue
                    if sa2i.fuzzySaOverlap(sa1j, dedupOffset):
                        markedForRemoval[i] = True
                elif sa1j.isLowQualFuzzy() and not sa2i.isLowQualFuzzy():
                    if markedForRemoval[j]:
                        continue
                    if sa1j.fuzzySaOverlap(sa2i, dedupOffset):
                        markedForRemoval[j] = True
            elif sa2i == sa1j:
                if sa1i.isLowQualFuzzy() and not sa2j.isLowQualFuzzy():
                    if markedForRemoval[i]:
                        continue
                    if sa1i.fuzzySaOverlap(sa2j, dedupOffset):
                        markedForRemoval[i] = True
                elif sa2j.isLowQualFuzzy() and not sa1i.isLowQualFuzzy():
                    if markedForRemoval[j]:
                        continue
                    if sa2j.fuzzySaOverlap(sa1i, dedupOffset):
                        markedForRemoval[j] = True
    return [resultsWithDups[i] for i in range(len(resultsWithDups)) if not markedForRemoval[i]]


def dedup(resultsWithDups, defaultReadLength):
    markedForRemoval = [False for _ in range(len(resultsWithDups))]
    for i in range(len(resultsWithDups)):
        lineChunksI = resultsWithDups[i]
        chr1i = lineChunksI[0]
        chr2i = lineChunksI[3]
        pos1i = int(lineChunksI[1])
        pos2i = int(lineChunksI[4])
        sa1i = SuppAlignment(lineChunksI[16])
        sa2i = SuppAlignment(lineChunksI[17])
        for j in range(len(resultsWithDups)):
            if i == j:
                continue
            lineChunksJ = resultsWithDups[j]
            chr1j = lineChunksJ[0]
            chr2j = lineChunksJ[3]
            pos1j = int(lineChunksJ[1])
            pos2j = int(lineChunksJ[4])
            sa1j = SuppAlignment(lineChunksJ[16])
            sa2j = SuppAlignment(lineChunksJ[17])
            if chr2i == chr2j and abs(pos2i - pos2j) < 2 * defaultReadLength:
                absoluteSuperiorityScore1 = absoluteSuperiority(chr1i, chr1j, pos1i, pos1j, sa1i, sa1j, defaultReadLength)
                if absoluteSuperiorityScore1 == 1:
                    markedForRemoval[j] = True
                elif absoluteSuperiorityScore1 == 2:
                    markedForRemoval[i] = True
            if chr1i == chr1j and abs(pos1i - pos1j) < 2 * defaultReadLength:
                absoluteSuperiorityScore2 = absoluteSuperiority(chr2i, chr2j, pos2i, pos2j, sa2i, sa2j, defaultReadLength)
                if absoluteSuperiorityScore2 == 1:
                    markedForRemoval[j] = True
                elif absoluteSuperiorityScore2 == 2:
                    markedForRemoval[i] = True
            if chr2i == chr1j and abs(pos2i - pos1j) < 2 * defaultReadLength:
                absoluteSuperiorityScoreCross1 = absoluteSuperiority(chr1i, chr2j, pos1i, pos2j, sa1i, sa2j, defaultReadLength)
                if absoluteSuperiorityScoreCross1 == 1:
                    markedForRemoval[j] = True
                elif absoluteSuperiorityScoreCross1 == 2:
                    markedForRemoval[i] = True
            if chr1i == chr2j and abs(pos1i - pos2j) < 2 * defaultReadLength:
                absoluteSuperiorityScoreCross2 = absoluteSuperiority(chr2i, chr1j, pos2i, pos1j, sa2i, sa1j, defaultReadLength)
                if absoluteSuperiorityScoreCross2 == 1:
                    markedForRemoval[j] = True
                elif absoluteSuperiorityScoreCross2 == 2:
                    markedForRemoval[i] = True
    return [resultsWithDups[i] for i in range(len(resultsWithDups)) if not markedForRemoval[i]]


initialResults = []
with open(sophiaOutput) as f:
    print(next(f).rstrip())
    for line in f:
        lineChunks = line.rstrip().split('\t')
        if int(lineChunks[9]) < 3:
            continue
        initialResults.append(lineChunks)

dedupFuzzyResults = dedupFuzzy(initialResults, defaultReadLength * 2)
dedupResults = dedup(dedupFuzzyResults, defaultReadLength)
for chunks in dedupResults:
    print(*chunks, sep='\t')
