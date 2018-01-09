#!/usr/bin/env python
from sys import argv

rightSideFile=argv[1]

with open(rightSideFile) as f:
    for line in f:
        lineChunks=line.rstrip().split('\t')
        geneMatches=",".join(set(lineChunks[4].split(',')))
        lineIndices=set(lineChunks[3].split(','))
        for index in lineIndices:
            print("\t".join(lineChunks[0:3] + [index,geneMatches] + lineChunks[5:]))
