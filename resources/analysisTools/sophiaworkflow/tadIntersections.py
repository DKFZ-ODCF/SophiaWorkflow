#!/usr/bin/env python
import numpy as np
from sys import argv,stderr
from itertools import chain

inputFile=argv[1]
inputFileWhitelist=argv[2]
consensusTadFile=argv[3]
smallEventThreshold=argv[4]
chromosomeList=[str(x) for x in range(1,23)]+["X","Y"]

whitelistedLineIndices=set()
with open(inputFileWhitelist) as f:
    for line in f:
        whitelistedLineIndices.add(int(line.rstrip()))

chrMidlines={"1":125000000,"2":93300000,"3":91000000,"4":50400000,"5":48400000,
             "6":61000000,"7":59900000,"8":45600000,"9":49000000,"10":40200000,
             "11":53700000,"12":35800000,"13":17900000,"14":17600000,"15":19000000,
             "16":36600000,"17":24000000,"18":17200000,"19":26500000,"20":27500000,
             "21":13200000,"22":14700000,"X":60600000,"Y":12500000}


class BpDigitizer:
    def __init__(self):
        self.windowStarts=dict()
        self.windowEnds=dict()
        self.genes=dict()
        self.cancerGenes=dict()
        self.lineIndices=dict()
        with open(consensusTadFile) as f:
            for line in f:
                lineChunks=line.rstrip().split('\t')
                chromosome=lineChunks[0]
                if chromosome not in self.windowStarts:
                    self.windowStarts[chromosome]=list()
                    self.windowEnds[chromosome]=list()
                    self.genes[chromosome]=list()
                    self.cancerGenes[chromosome]=list()
                    self.lineIndices[chromosome]=list()                    
                startPos=int(lineChunks[1])                
                endPos=int(lineChunks[2])
                self.windowStarts[chromosome].append(startPos)
                self.windowEnds[chromosome].append(endPos)
                self.genes[chromosome].append(lineChunks[3].split(','))
                self.cancerGenes[chromosome].append(lineChunks[4].split(','))
                self.lineIndices[chromosome].append(lineChunks[5])
        for chromosome in self.windowStarts:
            self.windowStarts[chromosome]=np.array(self.windowStarts[chromosome])
            self.windowEnds[chromosome]=np.array(self.windowEnds[chromosome])
    def processIntrachromosomalBps(self,bpChr,bpPos1,bpPos2,exclusionCandidacy,lowQual):
        eventSize=bpPos2-bpPos1
        windowOccupancy1=((bpPos1>=self.windowStarts[bpChr]) & (bpPos1<=self.windowEnds[bpChr]))
        windowOccupancy2=((bpPos2>=self.windowStarts[bpChr]) & (bpPos2<=self.windowEnds[bpChr]))
        windowOccupancy=(windowOccupancy1 | windowOccupancy2)
        iFirst=-1
        smallBorderHit=False
        if np.sum(windowOccupancy)==2:
            for i,x in enumerate(windowOccupancy):
                if x:
                    if iFirst==-1:
                        iFirst=i
                    else:
                        if i==iFirst+1:
                            smallBorderHit=True
                        break
        firstHit=0
        lastHit=0
        for x in range(len(windowOccupancy)):
            if windowOccupancy[x]:
                firstHit=x
                break
        for x in reversed(range(len(windowOccupancy))):
            if windowOccupancy[x]:
                lastHit=x
                break
        if eventSize < 10e6 or ((lastHit-firstHit)<=3 and not(bpPos1<=chrMidlines[bpChr] and bpPos2>=chrMidlines[bpChr])):
            for x in range(firstHit,lastHit+1):
                windowOccupancy[x]=True
        initialHits={x for x in range(len(windowOccupancy)) if windowOccupancy[x]}
        if bpChr not in {"hs37d5", "MT"}:
            if (not smallBorderHit and np.sum(windowOccupancy)>1) or eventSize>50000:
                if 0 in initialHits:
                    if (not (self.windowEnds[bpChr][0]<=chrMidlines[bpChr] and self.windowStarts[bpChr][1]>=chrMidlines[bpChr])) and (abs(self.windowStarts[bpChr][1] - bpPos2) < 2e6):
                        windowOccupancy[1]=True
                if len(windowOccupancy)-1 in initialHits:
                    if not (self.windowEnds[bpChr][-2]<=chrMidlines[bpChr] and self.windowStarts[bpChr][-1]>=chrMidlines[bpChr]) and (abs(bpPos1 - self.windowEnds[bpChr][-2]) < 2e6):
                        windowOccupancy[-2]=True
                for i in range(1,len(windowOccupancy)-1):
                    if i in initialHits:
                        if not (self.windowEnds[bpChr][i-1]<=chrMidlines[bpChr] and self.windowStarts[bpChr][i]>=chrMidlines[bpChr]) and (abs(self.windowStarts[bpChr][i] - bpPos2) < 2e6):
                            windowOccupancy[i-1]=True
                        if not (self.windowEnds[bpChr][i]<=chrMidlines[bpChr] and self.windowStarts[bpChr][i+1]>=chrMidlines[bpChr]) and (abs(bpPos1-self.windowEnds[bpChr][i-1]) < 2e6):
                            windowOccupancy[i+1]=True
        indices=[self.lineIndices[bpChr][i] for i,x in enumerate(windowOccupancy) if x]
        indicesStr=""
        if not (lowQual or (np.sum(windowOccupancy)<2 and exclusionCandidacy)):
            if exclusionCandidacy and bpPos2-bpPos1<1000 and smallBorderHit:
                indicesStr=','.join([x+"*" for x in indices])
            else:
                indicesStr=','.join(indices)
        else:
            indicesStr=','.join([x+"*" for x in indices])
        tmpGenes=set(chain(*[self.genes[bpChr][i] for i,x in enumerate(windowOccupancy) if x]))
        if len(tmpGenes)>1 and "." in tmpGenes:
            tmpGenes.remove(".")
        genes=','.join(tmpGenes)
        tmpCancerGenes=set(chain(*[self.cancerGenes[bpChr][i] for i,x in enumerate(windowOccupancy) if x]))
        if len(tmpCancerGenes)>1 and "." in tmpCancerGenes:
            tmpCancerGenes.remove(".")
        cancerGenes=','.join(tmpCancerGenes)
        return [indicesStr,genes,cancerGenes]
    def processBp(self,bpChr,bpPos):
        windowOccupancy=((bpPos>=self.windowStarts[bpChr]) & (bpPos<=self.windowEnds[bpChr]))
        initialHits=set()
        for i in range(len(windowOccupancy)):
            if windowOccupancy[i]:
                initialHits.add(i)
        if bpChr not in {"hs37d5", "MT"}:
            if 0 in initialHits:
                if (not (self.windowEnds[bpChr][0]<=chrMidlines[bpChr] and self.windowStarts[bpChr][1]>=chrMidlines[bpChr])) and (abs(self.windowStarts[bpChr][1] - bpPos) < 2e6):
                    windowOccupancy[1]=True
            if len(windowOccupancy)-1 in initialHits:
                if (not (self.windowEnds[bpChr][-2]<=chrMidlines[bpChr] and self.windowStarts[bpChr][-1]>=chrMidlines[bpChr])) and (abs(bpPos-self.windowEnds[bpChr][-2]) < 2e6):
                    windowOccupancy[-2]=True
            for i in range(1,len(windowOccupancy)-1):
                if i in initialHits:
                    if (not (self.windowEnds[bpChr][i-1]<=chrMidlines[bpChr] and self.windowStarts[bpChr][i]>=chrMidlines[bpChr])) and (abs(self.windowStarts[bpChr][i] - bpPos) < 2e6):
                        windowOccupancy[i-1]=True
                    if (not (self.windowEnds[bpChr][i-1]<=chrMidlines[bpChr] and self.windowStarts[bpChr][i]>=chrMidlines[bpChr])) and (abs(bpPos-self.windowEnds[bpChr][i-1]) < 2e6):
                        windowOccupancy[i+1]=True
        indices=[self.lineIndices[bpChr][i] for i,x in enumerate(windowOccupancy) if x]
        genes=[self.genes[bpChr][i] for i,x in enumerate(windowOccupancy) if x]
        cancerGenes=[self.cancerGenes[bpChr][i] for i,x in enumerate(windowOccupancy) if x]
        return [indices,genes,cancerGenes]
    def processResults(self,svResults):
        with open(svResults) as f:
            for lineIndex,svLine in enumerate(f):
                lineChunks=svLine.rstrip().split('\t')
                eventScore=int(lineChunks[9])
                lowQualCandidate=eventScore<3
                chr1=lineChunks[0]
                chr2=lineChunks[3]
                coord1=int(lineChunks[1])
                coord2=int(lineChunks[4])
                exclusionCandidate=False
                if chr1 == chr2:
                    if chr1=="hs37d5":
                        exclusionCandidate=True
                    else:
                        if sorted(lineChunks[18])==sorted(lineChunks[28]):
                            eventSize= abs(coord1-coord2)
                            if eventSize<5000:
                                if lineChunks[18]==".":
                                    if not (lineIndex in whitelistedLineIndices):
                                        exclusionCandidate=True
                                else:
                                    if "intron" in lineChunks[18] and not lineChunks[18].endswith("_intron1"):
                                        if not (lineIndex in whitelistedLineIndices):
                                            exclusionCandidate=True
                    minPos=min(int(lineChunks[1]),int(lineChunks[4]))
                    maxPos=max(int(lineChunks[1]),int(lineChunks[4]))
                    indicesStr,genes,cancerGenes=self.processIntrachromosomalBps(lineChunks[0], minPos, maxPos,exclusionCandidate,lowQualCandidate)
                    print("TRUE" if (lineIndex in whitelistedLineIndices and maxPos-minPos<5000) else "FALSE",indicesStr,genes,cancerGenes,sep='\t')
                else:
                    res=[[],[],[]]
                    item1,item2,item3=self.processBp(lineChunks[0], int(lineChunks[1]))
                    res[0]+=item1
                    res[1]+=item2
                    res[2]+=item3
                    item1,item2,item3=self.processBp(lineChunks[3], int(lineChunks[4]))
                    res[0]+=item1
                    res[1]+=item2
                    res[2]+=item3
                    indicesStr=""
                    if lowQualCandidate:
                        indicesStr=','.join([x+"*" for x in res[0]])
                    else:
                        indicesStr=','.join(res[0])
                    tmpGenes=set(chain(*res[1]))
                    if len(tmpGenes)>1 and "." in tmpGenes:
                        tmpGenes.remove(".")
                    genes=','.join(tmpGenes)
                    tmpCancerGenes=set(chain(*res[2]))
                    if len(tmpCancerGenes)>1 and "." in tmpCancerGenes:
                        tmpCancerGenes.remove(".")
                    cancerGenes=','.join(tmpCancerGenes)
                    print("FALSE",indicesStr,genes,cancerGenes,sep='\t')
tmpObj=BpDigitizer()
tmpObj.processResults(inputFile)
