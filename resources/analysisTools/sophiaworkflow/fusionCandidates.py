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
# Add the columns "directFusionCandidates", "directFusionCandidatesBothCancer",
# "indirectFusionCandidatesLeftCancerRightAny", "indirectFusionCandidatesRightCancerLeftAny",
# and "indirectFusionCandidatesAny" to the input TSV.
#
# Note that some output fields only contain the first gene name from a list in one of the input colums. If these
# input lists have quasi random order (e.g. being just the random-ordered fields of a dictionary) then these
# columns may appear to contain unpredictable values -- in particular different values for different runs
# of the workflow.
#
import fileinput
import re
import itertools
import sys


MAXDISTANCE = 2000000


def natural_string_sort(strings: list) -> list:
    """Sort a list of strings lowercase-alphabetically, for the alphabetical character parts, and numerically, for the
    numeric parts."""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(strings, key=alphanum_key)


def firstRawGeneName(extendedGeneName: str) -> str:
    """Extract rawName from in input string of the form 'rawName(|restA)?(;restB)*'. The ;-separated expression will
    probably be also extended gene names, but this function only returns the *first* gene name it can parse from the
    input string. If 'rawName(|restA)' contains an underscore, the last underscore-separated component will be dropped
    and 'rawName' (left of leftmost '|') will be returned. Thus if the rightmost '_' is left of the '|' the rawName
    will be shorter than the part left of the '|'."""
    extendedGeneName = extendedGeneName.split(';')[0]
    if '_' in extendedGeneName:
        return '_'.join(extendedGeneName.split('_')[:-1]).split('|')[0]
    else:
        return extendedGeneName.split('|')[0]


def getFusionList(unfilteredGenes1: list, unfilteredGenes2: list) -> list:
    """Given two lists of genes, genes1 and genes2, return a list of gene pairs (lists) [g1, g2]. The g1 and g2 values
    will be raw gene names, unless the raw gene names g1 == g2, in which case the function falls back to returning
    extended gene names. Note that because extended gene names that are both in genes1 and genes2 are currently removed
    from the input, always g1 != g2 (raw or extended). Note that g1 will be derived from genes1 and g2 from genes2."""
    genes1 = [x for x in unfilteredGenes1 if x not in unfilteredGenes2]
    genes2 = [x for x in unfilteredGenes2 if x not in unfilteredGenes1]
    fusionList = list()
    for x in set([x for x in itertools.product(genes1, genes2)]):
        rawName1 = firstRawGeneName(x[0])
        rawName2 = firstRawGeneName(x[1])
        if rawName1 == rawName2:
            fusionList.append(list(x))
        else:
            fusionList.append([rawName1, rawName2])
    return fusionList


def deduplicateFusionsList(fusionList: list) -> list:
    """Identify pairs (a, b) and (b, a). Note that the input pairs (x, y) are returned as (x, y) and never (y, x),
    i.e. the order (x, y) is preserved. Furthermore, of a pair (a, b) and a pair (b, a) following later in the list,
    only the first (a, b) is returned, but (b, a) is dropped."""
    deduplicatedFusions = list()
    seenFusions = set()
    for fusion in fusionList:
        # FusionIDs are (naturally) sorted to yield unique identifiers and remove symmetric pairs
        # (a,b)/(b,a) from the fusionList. Note that fusionID is thrown away in the end, because only the original
        # `fusion` variable is added to the output list. Thus, the sorting is only for the deduplication.
        fusionID = '-'.join(natural_string_sort(fusion))
        if fusionID not in seenFusions:
            seenFusions.add(fusionID)
            deduplicatedFusions.append(fusion)
    return deduplicatedFusions


def getDirectFusionStr(gene1Raw: str, gene2Raw: str) -> str:
    """Given two strings of comma-separated extended gene names, produce a string representing direct fusions. The
    string consists of a comma-separated list of pairs g1-g2, where g1 derives from input gene1Raw and g2 derives from
    gene2Raw. g1 and g2 are usually raw gene names, unless (1) they occur in both lists (see getFusionList()), or (2)
    the corresponding raw gene names in the pair g1-g2 are identical, in which case a fallback to the extended gene name
    happens. Note that the output order is stable by simply sorting the pair strings "g1-g2" before output. If a gene is
    paired with '.', then, instead of '.', '(TRUNC)' is returned. If no pairs are found '.' is returned. This is for
    instance the case if gene1Raw or gene2Raw input is empty (i.e. '.' or '')."""
    if gene1Raw in {".", ""} and gene2Raw in {".", ""}:
        return "."
    if gene1Raw == gene2Raw:
        return "."
    genes1 = natural_string_sort(gene1Raw.split(','))
    genes2 = natural_string_sort(gene2Raw.split(','))
    fusionList = getFusionList(genes1, genes2)
    if len(fusionList) > 0:
        deduplicatedFusionList = deduplicateFusionsList(fusionList)
        fusionCandidates = list()
        for fusion in deduplicatedFusionList:
            if fusion[0] in {'', '.'}:
                fusion[0] = "(TRUNC)"
            if fusion[1] in {'', '.'}:
                fusion[1] = "(TRUNC)"
            if fusion[0] != fusion[1]:
                # Eventually, join to fused pairs to "g1-g2" pairs, with g1<-genes1 & g2<-genes2
                fusionCandidates.append('-'.join(fusion))
        if len(fusionCandidates) > 0:
            return ','.join(sorted(fusionCandidates))
        else:
            return "."
    else:
        return "."


def extractIdListFromListOfIdSemicolonElems(string: str) -> str:
    """Given identifier lists of the form "a;a1,b;b1,c;c1", extract the list "a,b,c"."""
    return ','.join([elem.split(';')[0] for elem in string.split(',')])


def distanceFilter(value: str, distanceStr: str, max_distance: int = MAXDISTANCE) -> str:
    if value not in {'', '.'}:
        if int(distanceStr) > max_distance:
            return ""
        else:
            return value
    else:
        return ""


def getIndirectFusionComponentStr(geneRaw: str, geneNearestDownstream: str, geneNearestUpstream: str) -> str:
    componentStr = "None"
    if geneRaw in {"", "."}:
        if (geneNearestDownstream not in {"", "."}) or (geneNearestUpstream not in {"", "."}):
            if geneNearestUpstream not in {"", "."}:
                componentStr += "~" + geneNearestUpstream
            componentStr += "/"
            if geneNearestDownstream not in {"", "."}:
                componentStr += "~" + geneNearestDownstream
        else:
            componentStr = "(TRUNC)"
    else:
        componentStr = geneRaw
    return componentStr


def getIndirectFusionStr(gene1: str, gene1NearestUpstream: str, gene1NearestDownstream: str,
                         gene2: str, gene2NearestUpstream: str, gene2NearestDownstream: str) -> str:
    leftComponent =  getIndirectFusionComponentStr(gene1, gene1NearestUpstream, gene1NearestDownstream)
    rightComponent = getIndirectFusionComponentStr(gene2, gene2NearestUpstream, gene2NearestDownstream)
    if leftComponent != rightComponent:
        return leftComponent + "-" + rightComponent
    else:
        return "."


line_no = 0
try:
    for line in fileinput.input():
        line_no += 1
        if line[0] == '#':
            print(line.rstrip(), "directFusionCandidates", "directFusionCandidatesBothCancer",
                  "indirectFusionCandidatesLeftCancerRightAny",
                  "indirectFusionCandidatesRightCancerLeftAny", "indirectFusionCandidatesAny", sep='\t')
        else:
            lineChunks = line.rstrip().split('\t')

            gene1Raw =                     extractIdListFromListOfIdSemicolonElems(lineChunks[20])
            gene1RawCancer =               extractIdListFromListOfIdSemicolonElems(lineChunks[21])
            gene1NearestUpstreamRaw =      distanceFilter(lineChunks[22],
                                                          distanceStr=lineChunks[23])
            gene1NearestUpstreamCancer =   distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[24]),
                                                          distanceStr=lineChunks[25])
            gene1NearestDownstreamRaw =    distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[26]),
                                                          distanceStr=lineChunks[27])
            gene1NearestDownstreamCancer = distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[28]),
                                                          distanceStr=lineChunks[29])

            gene2Raw =                     extractIdListFromListOfIdSemicolonElems(lineChunks[30])
            gene2RawCancer =               extractIdListFromListOfIdSemicolonElems(lineChunks[31])
            gene2NearestUpstreamRaw =      distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[32]),
                                                          distanceStr=lineChunks[33])
            gene2NearestUpstreamCancer =   distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[34]),
                                                          distanceStr=lineChunks[35])
            gene2NearestDownstreamRaw =    distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[36]),
                                                          distanceStr=lineChunks[37])
            gene2NearestDownstreamCancer = distanceFilter(extractIdListFromListOfIdSemicolonElems(lineChunks[38]),
                                                          distanceStr=lineChunks[39])

            indirectFusionCandidatesLeftCancerRightAny = \
                getIndirectFusionStr(gene1RawCancer, gene1NearestUpstreamCancer, gene1NearestDownstreamCancer,
                                     gene2Raw, gene2NearestUpstreamRaw, gene2NearestDownstreamRaw)

            indirectFusionCandidatesRightCancerLeftAny = \
                getIndirectFusionStr(gene1Raw, gene1NearestUpstreamRaw, gene1NearestDownstreamRaw,
                                     gene2RawCancer, gene2NearestUpstreamCancer, gene2NearestDownstreamCancer)

            indirectFusionCandidatesAny = \
                getIndirectFusionStr(gene1Raw, gene1NearestUpstreamRaw, gene1NearestDownstreamRaw,
                                     gene2Raw, gene2NearestUpstreamRaw, gene2NearestDownstreamRaw)

            directFusionCandidates = \
                getDirectFusionStr(gene1Raw, gene2Raw)
            directFusionCandidatesBothCancer = \
                getDirectFusionStr(gene1RawCancer, gene2RawCancer)

            print(line.rstrip(),
                  directFusionCandidates, directFusionCandidatesBothCancer,
                  indirectFusionCandidatesLeftCancerRightAny, indirectFusionCandidatesRightCancerLeftAny,
                  indirectFusionCandidatesAny, sep='\t')

except Exception as e:
    print("Error while processing line {}: {}".format(line_no, e), file=sys.stderr)
    sys.exit(1)
