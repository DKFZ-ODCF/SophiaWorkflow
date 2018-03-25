#!/ibios/tbi_cluster/13.1/x86_64/bin/python3.4
from sys import argv,stderr
from sys import exit as sysexit
inputPath=argv[1]
eventDict=dict()
with open(inputPath) as inputHandle:
    for line in inputHandle:
        lineChunks=line.rstrip().split('\t')
        direction=int(lineChunks[3])
        lineIndex=int(lineChunks[4])
        gene=lineChunks[-1].split('_')[0]
        if lineIndex not in eventDict:
            eventDict[lineIndex]=[set(),set()]
        eventDict[lineIndex][direction].add(gene)
    
resultDict=dict()
for key in eventDict:
    event=eventDict[key]
    if len(event[0])>0 and len(event[1])>0:
        intersectionGenes=event[0].intersection(event[1])
        if len(intersectionGenes)>0:
            for gene in intersectionGenes:
                if gene not in resultDict:
                    resultDict[gene]=0
                resultDict[gene]+=1
for gene in resultDict:
    print(gene,resultDict[gene],sep='\t')
print("#done")
sysexit(0)
