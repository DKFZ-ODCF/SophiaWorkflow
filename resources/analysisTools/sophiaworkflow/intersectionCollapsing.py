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
# Aggregate all annotations (genes, cancer genes, superenhancers) for rows with the same ID column value
# (1-based: column 4). The input is the output from `bedtools intersect` in sophiaAnnotateAbridgedCaller.sh.
#
# Specifically, the following columns are used (as 1-based indices)
#
# 4 = ID
# 5 = hit class, 1 = gene, 2 = cancer gene, 3 = superenhancers
# 9 = annotation
#
# Lines (hits) with the same key are aggregated by combining into three sets the genes, cancer genes and superenhancers.
# NOTE: Lines with the same ID value are **only** combined, if they occur in a continuous sequence. If they are
# interrupted by lines with other ID values they are considered distinct, and result in distint output lines. There are
# NO checks, whether IDs occur again after an interruption.
#
# For each (continuous sequence of an) ID, the three sets of genes, cancer genes and superenhancers are printed as one
# TSV row with one column for each set.
#
# Usage:
#         intersectionCollapsing.py ${name}directHits[12]Pre
#
from sys import argv, stderr


def set_join(the_set: set, sep=",", default="."):
    if len(the_set) != 0:
        return sep.join(sorted(list(the_set)))
    else:
        return default


def print_hits(gene_hits: set, cancer_gene_hits: set, super_enhancer_hits: set):
    print(set_join(gene_hits),
          set_join(cancer_gene_hits),
          set_join(super_enhancer_hits),
          sep='\t')


if len(argv) != 2:
    print("Incorrect arguments. Usage: intersectionCollapsing.py ${name}directHits[12]Pre", file=stderr)
    exit(1)


inputFile = argv[1]


lastID = ""
geneHits = set()
cancerGeneHits = set()
superEnhancerHits = set()

with open(inputFile) as f:
    lineCount = 0
    for line in f:
        lineCount += 1
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
    if lineCount > 0:
        print_hits(geneHits, cancerGeneHits, superEnhancerHits)
