from sys import argv
inputFile=argv[1]

lastID=""
geneHits=set()
cancerGeneHits=set()
superEnhancerHits=set()
with open(inputFile) as f:
    for line in f:
        lineChunks=line.rstrip().split('\t')
        ID=lineChunks[3]
        if ID!=lastID:
            if lastID!="":
                hitStr="."
                cancerHitStr="."
                superEnhancerHitStr="."
                if len(geneHits)!=0:
                    hitStr=','.join(geneHits)
                if len(cancerGeneHits)!=0:
                    cancerHitStr=','.join(cancerGeneHits)
                if len(superEnhancerHits)!=0:
                    superEnhancerHitStr=','.join(superEnhancerHits)
                print(hitStr,cancerHitStr,superEnhancerHitStr,sep='\t')
            lastID=ID
            geneHits=set()
            cancerGeneHits=set()
            superEnhancerHits=set()
        if lineChunks[4]!=".":
            if lineChunks[4]=="1":
                geneHits.add(lineChunks[8])
            elif lineChunks[4]=="2":
                cancerGeneHits.add(lineChunks[8])
            elif lineChunks[4]=="3":
                superEnhancerHits.add(lineChunks[8])
    hitStr="."
    cancerHitStr="."
    superEnhancerHitStr="."
    if len(geneHits)!=0:
        hitStr=','.join(geneHits)
    if len(cancerGeneHits)!=0:
        cancerHitStr=','.join(cancerGeneHits)
    if len(superEnhancerHits)!=0:
        superEnhancerHitStr=','.join(superEnhancerHits)
    print(hitStr,cancerHitStr,superEnhancerHitStr,sep='\t')
