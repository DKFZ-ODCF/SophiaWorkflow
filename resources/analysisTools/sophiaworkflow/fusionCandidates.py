#!/ibios/tbi_cluster/13.1/x86_64/bin/python3.4
import fileinput
import re
import itertools
from sys import stderr

MAXDISTANCE=2000000

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
def rawGeneName(extendedGeneName):
    if '_' in extendedGeneName:
        return '_'.join(extendedGeneName.split('_')[:-1])
    else:
        return extendedGeneName

def getDirectFusion(gene1Raw,gene2Raw):
    if gene1Raw in {".",""} and gene2Raw in {".",""}:
        return "."
    if gene1Raw==gene2Raw:
        return "."
    gene1ListPre=gene1Raw.split(',')
    gene2ListPre=gene2Raw.split(',')
    genes1=natural_sort([x for x in gene1ListPre if x not in gene2ListPre])
    genes2=natural_sort([x for x in gene2ListPre if x not in gene1ListPre]) 
    fusionList=list()
    for x in set(itertools.product(genes1,genes2)):
        rawName1=rawGeneName(x[0])
        rawName2=rawGeneName(x[1])
        if rawName1==rawName2:
            fusionList.append(x)
        else:
            fusionList.append([rawName1,rawName2])
    sortedFusions=set()
    filteredFusionList=list()
    fusionCandidates=list()                
    if len(fusionList)>0:
        for fusion in fusionList:
            fusionID='-'.join(natural_sort(fusion))
            if fusionID not in sortedFusions:
                sortedFusions.add(fusionID)
                filteredFusionList.append(fusion)
        fusionCandidates=list()
        for fusionPre in filteredFusionList:
            fusion=list(fusionPre)
            if fusion[0] in {'','.'} :
                fusion[0]="(TRUNC)"
            if fusion[1] in {'','.'}:
                fusion[1]="(TRUNC)"
            if fusion[0]!=fusion[1]:
                fusionCandidates.append('-'.join(fusion))
        if len(fusionCandidates)>0:
            return ','.join(fusionCandidates)
        else:
            return "."
    else:
        return "."

for line in fileinput.input():
    if line[0]=='#':
        print(line.rstrip(),"directFusionCandidates","directFusionCandidatesBothCancer","indirectFusionCandidatesLeftCancerRightAny","indirectFusionCandidatesRightCancerLeftAny","indirectFusionCandidatesAny",sep='\t')
    else:
        directFusionCandidates=""
        directFusionCandidatesBothCancer=""
        indirectFusionCandidatesLeftCancerRightAny=""
        indirectFusionCandidatesRightCancerLeftAny=""
        indirectFusionCandidatesAny=""
        lineChunks=line.rstrip().split('\t')
        gene1Raw=lineChunks[18]
        gene1RawCancer=lineChunks[19]
        gene1NearestUpstreamRaw=lineChunks[20]
        if gene1NearestUpstreamRaw not in {'','.'}:
            gene1NearestUpstreamDistance=int(lineChunks[21])
            if gene1NearestUpstreamDistance>MAXDISTANCE:
                gene1NearestUpstreamRaw=""
        else:
            gene1NearestUpstreamRaw=""
        gene1NearestUpstreamCancer=lineChunks[22]
        if gene1NearestUpstreamCancer not in {'','.'}:
            gene1NearestUpstreamCancerDistance=int(lineChunks[23])
            if gene1NearestUpstreamCancerDistance>MAXDISTANCE:
                gene1NearestUpstreamCancer=""
        else:
            gene1NearestUpstreamCancer=""
        gene1NearestDownstreamRaw=lineChunks[24]
        if gene1NearestDownstreamRaw not in {'','.'}:
            gene1NearestDownstreamDistance=int(lineChunks[25])
            if gene1NearestDownstreamDistance>MAXDISTANCE:
                gene1NearestDownstreamRaw=""
        else:
            gene1NearestDownstreamRaw=""
        gene1NearestDownstreamCancer=lineChunks[26]
        if gene1NearestDownstreamCancer not in {'','.'}:
            gene1NearestDownstreamCancerDistance=int(lineChunks[27])
            if gene1NearestDownstreamCancerDistance>MAXDISTANCE:
                gene1NearestDownstreamCancer=""
        else:
            gene1NearestDownstreamCancer=""
        gene2Raw=lineChunks[28]
        gene2RawCancer=lineChunks[29]
        gene2NearestUpstreamRaw=lineChunks[30]
        if gene2NearestUpstreamRaw not in {'','.'}:
            gene2NearestUpstreamDistance=int(lineChunks[31])
            if gene2NearestUpstreamDistance>MAXDISTANCE:
                gene2NearestUpstreamRaw=""
        else:
            gene2NearestUpstreamRaw=""
        gene2NearestUpstreamCancer=lineChunks[32]
        if gene2NearestUpstreamCancer not in {'','.'}:
            gene2NearestUpstreamCancerDistance=int(lineChunks[33])
            if gene2NearestUpstreamCancerDistance>MAXDISTANCE:
                gene2NearestUpstreamCancer=""
        else:
            gene2NearestUpstreamCancer=""
        gene2NearestDownstreamRaw=lineChunks[34]
        if gene2NearestDownstreamRaw not in {'','.'}:
            gene2NearestDownstreamDistance=int(lineChunks[35])
            if gene2NearestDownstreamDistance>MAXDISTANCE:
                gene2NearestDownstreamRaw=""
        else:
            gene2NearestDownstreamRaw=""
        gene2NearestDownstreamCancer=lineChunks[36]
        if gene2NearestDownstreamCancer not in {'','.'}:
            gene2NearestDownstreamCancerDistance=int(lineChunks[37])
            if gene2NearestDownstreamCancerDistance>MAXDISTANCE:
                gene2NearestDownstreamCancer=""
        else:
            gene2NearestDownstreamCancer=""
        directFusionCandidates=getDirectFusion(gene1Raw,gene2Raw)
        directFusionCandidatesBothCancer=getDirectFusion(gene1RawCancer,gene2RawCancer)
        leftComponent=""
        if gene1RawCancer in {"","."}:
            if (gene1NearestDownstreamCancer not in {"","."}) or (gene1NearestUpstreamCancer not in {"","."}):
                if gene1NearestUpstreamCancer not in {"","."}:
                    leftComponent+="~"+gene1NearestUpstreamCancer
                leftComponent+="/"
                if gene1NearestDownstreamCancer not in {"","."}:
                    leftComponent+="~"+gene1NearestDownstreamCancer
            else:
                leftComponent="(TRUNC)"
        else:
            leftComponent=gene1RawCancer
        rightComponent=""
        if gene2Raw in {"","."}:
            if (gene2NearestUpstreamRaw not in {"","."}) or (gene2NearestDownstreamRaw not in {"","."}):
                if gene2NearestUpstreamRaw not in {"","."}:
                    rightComponent+="~"+gene2NearestUpstreamRaw
                rightComponent+="/"
                if gene2NearestDownstreamRaw not in {"","."}:
                    rightComponent+="~"+gene2NearestDownstreamRaw
            else:
                rightComponent="(TRUNC)"
        else:
            rightComponent=rawGeneName(gene2Raw)
        if leftComponent!=rightComponent:
            indirectFusionCandidatesLeftCancerRightAny=leftComponent+"-"+rightComponent
        else:
            indirectFusionCandidatesLeftCancerRightAny="."
        leftComponent=""
        if gene1Raw in {"","."}:
            if (gene1NearestUpstreamRaw not in {"","."}) or (gene1NearestDownstreamRaw not in {"","."}):
                if gene1NearestUpstreamRaw not in {"","."}:
                    leftComponent+="~"+gene1NearestUpstreamRaw
                leftComponent+="/"
                if gene1NearestDownstreamRaw not in {"","."}:
                    leftComponent+="~"+gene1NearestDownstreamRaw
            else:
                leftComponent="(TRUNC)"
        else:
            leftComponent=rawGeneName(gene1Raw)
        rightComponent=""
        if gene2RawCancer in {"","."}:
            if (gene2NearestUpstreamCancer not in {"","."}) or (gene2NearestDownstreamCancer not in {"","."}):
                if gene2NearestUpstreamCancer not in {"","."}:
                    rightComponent+="~"+gene2NearestUpstreamCancer
                rightComponent+="/"
                if gene2NearestDownstreamCancer not in {"","."}:
                    rightComponent+="~"+gene2NearestDownstreamCancer
            else:
                rightComponent="(TRUNC)"
        else:
            rightComponent=gene2RawCancer
        if leftComponent!=rightComponent:
            indirectFusionCandidatesRightCancerLeftAny=leftComponent+"-"+rightComponent
        else:
            indirectFusionCandidatesRightCancerLeftAny="."
        leftComponent=""
        
        if gene1Raw in {"","."}:
            if (gene1NearestUpstreamRaw not in {"","."}) or (gene1NearestDownstreamRaw not in {"","."}):
                if gene1NearestUpstreamRaw not in {"","."}:
                    leftComponent+="~"+gene1NearestUpstreamRaw
                leftComponent+="/"
                if gene1NearestDownstreamRaw not in {"","."}:
                    leftComponent+="~"+gene1NearestDownstreamRaw
            else:
                leftComponent="(TRUNC)"
        else:
            leftComponent=rawGeneName(gene1Raw)
        rightComponent=""
        if gene2Raw in {"","."}:
            if (gene2NearestUpstreamRaw not in {"","."}) or (gene2NearestDownstreamRaw not in {"","."}):
                if gene2NearestUpstreamRaw not in {"","."}:
                    rightComponent+="~"+gene2NearestUpstreamRaw
                rightComponent+="/"
                if gene2NearestDownstreamRaw not in {"","."}:
                    rightComponent+="~"+gene2NearestDownstreamRaw
            else:
                rightComponent="(TRUNC)"
        else:
            rightComponent=rawGeneName(gene2Raw)
        if leftComponent!=rightComponent:
            indirectFusionCandidatesAny=leftComponent+"-"+rightComponent
        else:
            indirectFusionCandidatesAny="."        
        print(line.rstrip(),directFusionCandidates,directFusionCandidatesBothCancer,indirectFusionCandidatesLeftCancerRightAny,indirectFusionCandidatesRightCancerLeftAny, indirectFusionCandidatesAny,sep='\t')
