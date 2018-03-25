#!/ibios/tbi_cluster/13.1/x86_64/bin/python3.4
from sys import argv

righSideFile=argv[1]

with open(righSideFile) as f:
    for line in f:
        lineChunks=line.rstrip().split('\t')
        geneMatches=",".join(set(lineChunks[4].split(',')))
        lineIndices=set(lineChunks[3].split(','))
        for index in lineIndices:
            print("\t".join(lineChunks[0:3] + [index,geneMatches] + lineChunks[5:]))
