#!/usr/bin/env python
from sys import argv,stderr
inputPath=argv[1]
lineIndex=0
with open(inputPath) as inputHandle:
    for line in inputHandle:
        if line[0]!='#':
            lineChunks=line.rstrip().split('\t')
            chr1=lineChunks[0]
            chr2=lineChunks[3]
            eventType=lineChunks[8]
            lineChunks[4]=int(lineChunks[4])+1
            lineChunks[5]=int(lineChunks[5])+1
            genes1=list()
            genes2=list()
            eventScore=int(lineChunks[9])
            if eventType=="DEL" or (eventType=="TRA" and chr1==chr2):
                print(lineChunks[0],lineChunks[1],lineChunks[2],0,lineIndex,sep='\t')
                print(lineChunks[3],lineChunks[4],lineChunks[5],1,lineIndex,sep='\t')
            lineIndex+=1
