from sys import argv
inputBedpe=argv[1]
with open(inputBedpe) as f:
    for i,line in enumerate(f):
        lineChunks=line.rstrip().split('\t')
        if lineChunks[0]!=lineChunks[3]:
            print(lineChunks[0],int(9e9),int(9e9+1),i,sep='\t')
        else:
            minPos,maxPos=sorted([int(lineChunks[1]),int(lineChunks[4])])
            print(lineChunks[0],minPos,maxPos,i,sep='\t')