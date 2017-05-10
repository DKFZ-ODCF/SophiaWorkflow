from sys import argv,stderr
from itertools import product
preFilteredBedpe=argv[1]

broadProd=list(product(range(51),range(51)))
narrowProd=list(product(range(10),range(10)))

lineIndex=0
with open(preFilteredBedpe) as inputHandle:
    for line in inputHandle:
        if line[0]!='#':
            lineChunks=line.rstrip().split('\t')
            eventType=lineChunks[8]
            if lineChunks[0]==lineChunks[3] and eventType in {"TRA","DEL"}:
                minPos=min(int(lineChunks[1]),int(lineChunks[4]))
                maxPos=max(int(lineChunks[1]),int(lineChunks[4]))
                eventScore=int(lineChunks[9])
                for iterPair in broadProd:
                    print(lineChunks[0],minPos-iterPair[0],maxPos-1+iterPair[1],lineIndex,sep='\t')
            lineIndex+=1
