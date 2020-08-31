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
#
#
# Aggregate all annotations for rows with the same ID column value (1-based: column 4).
#
# Usage:
#         intersectionCollapsing.py ${name}directHits[12]Pre
#
# The input file comes from bedtools intersect. See sophiaAnnotateAbridgedCaller.sh for details..
#
from sys import argv


def set_join(the_set: set, sep=",", default="."):
    if len(the_set) != 0:
        return sep.join(the_set)
    else:
        return default


def print_hits(gene_hits: set, cancer_gene_hits: set, super_enhancer_hits: set):
    print(set_join(gene_hits),
          set_join(cancer_gene_hits),
          set_join(super_enhancer_hits),
          sep='\t')


inputFile = argv[1]


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
                print_hits(geneHits, cancerGeneHits, superEnhancerHits)

            lastID = ID
            geneHits = set()
            cancerGeneHits = set()
            superEnhancerHits = set()

        hitClass = lineChunks[4]
        if hitClass != ".":
            annotation = lineChunks[8]
            if hitClass == "1":
                geneHits.add(annotation)
            elif hitClass == "2":
                cancerGeneHits.add(annotation)
            elif hitClass == "3":
                superEnhancerHits.add(annotation)
            else:
                raise "Unknown hit class (.=unclassified, 1=gene, 2=cancerGene, 3=superEnhancer): '%'".\
                    format(hitClass)

    print_hits(geneHits, cancerGeneHits, superEnhancerHits)
