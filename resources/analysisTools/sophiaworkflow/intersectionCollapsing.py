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

inputFile = argv[1]


def print_line(hit_str, cancer_hit_str, super_enhancer_hit_str):
    """Print a line of information unless it is trivial/all '.'"""
    if hit_str != "." or cancer_hit_str != "." or super_enhancer_hit_str != ".":
        print(hit_str, cancer_hit_str, super_enhancer_hit_str, sep='\t')


lastID = ""
geneHits = set()
cancerGeneHits = set()
superEnhancerHits = set()
with open(inputFile) as f:
    for line in f:
        lineChunks = line.rstrip().split('\t')
        ID = lineChunks[3]
        if ID != lastID:
            if lastID != "":
                hitStr = "."
                cancerHitStr = "."
                superEnhancerHitStr = "."
                if len(geneHits) != 0:
                    hitStr = ','.join(geneHits)
                if len(cancerGeneHits) != 0:
                    cancerHitStr = ','.join(cancerGeneHits)
                if len(superEnhancerHits) != 0:
                    superEnhancerHitStr = ','.join(superEnhancerHits)
                print_line(hitStr, cancerHitStr, superEnhancerHits)
            lastID = ID
            geneHits = set()
            cancerGeneHits = set()
            superEnhancerHits = set()
        if lineChunks[4] != ".":
            if lineChunks[4] == "1":
                geneHits.add(lineChunks[8])
            elif lineChunks[4] == "2":
                cancerGeneHits.add(lineChunks[8])
            elif lineChunks[4] == "3":
                superEnhancerHits.add(lineChunks[8])
    hitStr = "."
    cancerHitStr = "."
    superEnhancerHitStr = "."
    if len(geneHits) != 0:
        hitStr = ','.join(geneHits)
    if len(cancerGeneHits) != 0:
        cancerHitStr = ','.join(cancerGeneHits)
    if len(superEnhancerHits) != 0:
        superEnhancerHitStr = ','.join(superEnhancerHits)
    print_line(hitStr, cancerHitStr, superEnhancerHitStr)
