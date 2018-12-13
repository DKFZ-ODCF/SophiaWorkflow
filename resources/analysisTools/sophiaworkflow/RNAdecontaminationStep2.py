#!/usr/bin/env python
from sys import argv,stderr
from itertools import product
preFilteredBedpe=argv[1]
lineIndicesFile=argv[2]

def getIntronExonIndex(geneRawName):
    geneChunks=geneRawName.split('_')
    geneName="."
    if geneChunks[0]!=".":
        geneName='_'.join(geneChunks[:-1])
    index=1
    intronStatus="intron" in geneChunks[-1] or geneName=="."
    if geneChunks[-1]!=".":
        index=int(geneChunks[-1].split('intron')[-1])
    return [geneName,index,intronStatus]

lineIndices=set()
with open(lineIndicesFile) as f:
    for line in f:
        lineIndex=int(line.rstrip())
        lineIndices.add(lineIndex)

lineIndex=-1
with open(preFilteredBedpe) as inputHandle:
    for line in inputHandle:
        if line[0]!='#':
            lineIndex+=1
            skipLine=False
            if lineIndex not in lineIndices:
                lineChunks=line.rstrip().split('\t')
                eventType=lineChunks[8]
                eventScore=int(lineChunks[9])
                if lineChunks[0]==lineChunks[3] and lineChunks[11]!="INV":
                    gene1Raw=lineChunks[18]
                    gene2Raw=lineChunks[28]
                    gene1ListPre=[getIntronExonIndex(x) for x in gene1Raw.split(',') if '_' in x]
                    gene2ListPre=[getIntronExonIndex(x) for x in gene2Raw.split(',') if '_' in x]
                    for eventPair in product(gene1ListPre,gene2ListPre):
                        if eventPair[0][0]==eventPair[1][0]:
                            if eventPair[0][0]!=".":
                                if (eventPair[0][2] or eventPair[1][2]):
                                    if abs(eventPair[0][1]-eventPair[1][1])<2:
                                        skipLine=True
                                        break
                                else:
                                    if eventType=="DEL" and abs(eventPair[0][1]-eventPair[1][1])==1:
                                        skipLine=True
                                        break
                if not skipLine:
                    print(line.rstrip())
        else:
            print(line.rstrip())
