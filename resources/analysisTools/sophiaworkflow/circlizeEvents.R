#!/usr/bin/env Rscript

library(circlize)

textPdf <- function(filename, text) {
    pdf(filename)
    plot.new()
    mtext(text)
    dev.off()
}

args=commandArgs(TRUE)
annosFile=args[1]
pid=args[2]
scoreThreshold=as.numeric(args[3])

annos=read.delim(annosFile,header = TRUE,stringsAsFactors = FALSE)

scaledOutput=paste(annosFile,"score",scoreThreshold,"scaled.pdf",sep="_")
if(nrow(annos) == 0){
    textPdf(scaledOutput, paste(basename(scaledOutput), "No data to plot.", sep=":\n"))
    stop()
}

annos[,1]=as.character(annos[,1])
annos[,4]=as.character(annos[,4])
annos$chrom1PreDecoyRemap=as.character(annos$chrom1PreDecoyRemap)
annos$chrom2PreDecoyRemap=as.character(annos$chrom2PreDecoyRemap)

annos=annos[annos$eventScore>=scoreThreshold,]
annos=annos[!(annos$chrom1PreDecoyRemap=="hs37d5" & annos$chrom2PreDecoyRemap=="hs37d5"),]
annos=annos[!(annos$chrom1=="hs37d5" & annos$chrom2PreDecoyRemap=="hs37d5"),]
annos=annos[!(annos$chrom1PreDecoyRemap=="hs37d5" & annos$chrom2=="hs37d5"),]
annos=annos[!(annos[,1]=="X" & annos[,4]=="Y"),]
annos=annos[!(annos[,1]=="Y" & annos[,4]=="X"),]

annos=annos[!((startsWith(annos[,1],"h") |startsWith(annos[,1],"G") | startsWith(annos[,1],"N") | startsWith(annos[,1],"H"))  & (startsWith(annos[,4],"h") | startsWith(annos[,4],"G") | startsWith(annos[,4],"N") | startsWith(annos[,4],"H"))),]

annos=annos[!(annos[,1]=="hs37d5" & grepl("UNKNOWN",annos[,8])),]
annos=annos[!(startsWith(annos[,1],"G") & grepl("UNKNOWN",annos[,8])),]
annos=annos[!(startsWith(annos[,1],"N") & grepl("UNKNOWN",annos[,8])),]

if(nrow(annos) == 0){
    textPdf(scaledOutput, paste(basename(scaledOutput), "No data to plot.", sep=":\n"))
    stop()
}

annos$dummyChromosome=rep(FALSE,nrow(annos))

mtLocs=annos[,4]=="MT"
if(any(mtLocs)){
  annos[mtLocs,4]=annos[mtLocs,1]
  annos[mtLocs,5]=annos[mtLocs,3]
  annos[mtLocs,6]=annos[mtLocs,3]+1
  annos$dummyChromosome[mtLocs]=TRUE  
}

decoyLocs=annos[,4]=="hs37d5" | startsWith(annos[,4],"G") | startsWith(annos[,4],"N")
if(any(decoyLocs)){
  annos[decoyLocs,4]=annos[decoyLocs,1]
  annos[decoyLocs,5]=annos[decoyLocs,3]
  annos[decoyLocs,6]=annos[decoyLocs,3]+1
  annos$svtype[decoyLocs]="UNKNOWN_DECOY"
  annos$dummyChromosome[decoyLocs]=TRUE
}
decoyLocs2=annos[,1]=="hs37d5" | startsWith(annos[,1],"G") | startsWith(annos[,1],"N")
if(any(decoyLocs2)){
	annos[decoyLocs2,1]=annos[decoyLocs2,4]
	annos[decoyLocs2,2]=annos[decoyLocs2,6]
	annos[decoyLocs2,3]=annos[decoyLocs2,6]+1
	annos$svtype[decoyLocs2]="UNKNOWN_DECOY"
	annos$dummyChromosome[decoyLocs2]=TRUE
}

telInsLocs=(annos[,4]=="5" & annos[,5] < 11890)
if(any(telInsLocs)){
	annos[telInsLocs,4]=annos[telInsLocs,1]
	annos[telInsLocs,5]=annos[telInsLocs,2]
	annos[telInsLocs,6]=annos[telInsLocs,3]+1
	annos$svtype[telInsLocs]="TELINS"
	annos$dummyChromosome[telInsLocs]=TRUE
}


annos$clonality=rep(0,nrow(annos))
unknownLocs1=annos$clonalityRatio1=="UNKNOWN"
annos$clonalityRatio1[unknownLocs1]=0.01
annos$clonalityRatio1=as.numeric(annos$clonalityRatio1)
annos$clonalityRatio2=as.character(annos$clonalityRatio2)
unknownLocs=annos$clonalityRatio2=="UNKNOWN"
annos$clonalityRatio2[unknownLocs]=NA
annos$clonalityRatio2=as.numeric(annos$clonalityRatio2)
imbalancedLocs=!unknownLocs&abs(annos$clonalityRatio1-annos$clonalityRatio2)>0.3
balancedLocs=!unknownLocs&!imbalancedLocs
if(any(unknownLocs)){
  annos$clonality[unknownLocs]=annos$clonalityRatio1[unknownLocs]
}
if(any(imbalancedLocs)){
  annos$clonality[imbalancedLocs]=pmax(annos$clonalityRatio1[imbalancedLocs],annos$clonalityRatio2[imbalancedLocs])
}
if(any(balancedLocs)){
  annos$clonality[balancedLocs]=0.5*(annos$clonalityRatio1[balancedLocs]+annos$clonalityRatio2[balancedLocs])
}

annos$scaledClonality=annos$clonality/(quantile((annos$clonality))[[4]])
annos$scaledClonality[annos$scaledClonality > 1] <- 1

coords1=annos[,1:3]
coords1[,2]=as.numeric(coords1[,2])
coords1[,1]=paste0("chr",coords1[,1])
coords2=annos[,4:6]
coords2[,2]=as.numeric(coords2[,2])
coords2[,1]=paste0("chr",coords2[,1])

colnames(coords1)=c("chr","start","end")
colnames(coords2)=c("chr","start","end")
eventTypes=annos$svtype
eventInversion=annos$eventInversion
newInversions=eventTypes!="UNKNOWN_DECOY"&eventTypes!="TELINS"&eventTypes!="TRA" & eventTypes!="INV" & eventInversion=="INV"
eventTypes[newInversions]="INV"

smallEventSelection=!annos$dummyChromosome & !is.na(annos$eventSize) & annos$eventSize < 9e6
mediumEventSelection=!annos$dummyChromosome & !is.na(annos$eventSize) & annos$eventSize >= 9e6 & annos$eventSize < 18e6
largeEventSelection=is.na(annos$eventSize) | annos$dummyChromosome | annos$eventSize >= 18e6

eventColors=rep("gray",length(eventTypes))
eventScaledColors=rep("gray",length(eventTypes))
smallEventScaling=ifelse(smallEventSelection,0.5,1)
alphas=annos$clonality*smallEventScaling
scaledAlphas=annos$scaledClonality*smallEventScaling

invLocs=eventTypes=="INV"
delLocs=eventTypes=="DEL"
dupLocs=eventTypes=="DUP"
traLocs=eventTypes=="TRA"
telInsLocs=eventTypes=="TELINS"
otherLocs=eventTypes=="UNKNOWN" | eventTypes=="UNKNOWN_DECOY"

if(any(invLocs)){
eventColors[invLocs]=rgb(t(col2rgb("black"))/255, alpha = alphas[invLocs])
eventScaledColors[invLocs]=rgb(t(col2rgb("black"))/255, alpha = scaledAlphas[invLocs])  
}
if(any(delLocs)){
eventColors[delLocs]=rgb(t(col2rgb("blue"))/255, alpha = alphas[delLocs])
eventScaledColors[delLocs]=rgb(t(col2rgb("blue"))/255, alpha = scaledAlphas[delLocs])  
}
if(any(dupLocs)){
eventColors[dupLocs]=rgb(t(col2rgb("red"))/255, alpha = alphas[dupLocs])
eventScaledColors[dupLocs]=rgb(t(col2rgb("red"))/255, alpha = scaledAlphas[dupLocs])  
}
if(any(traLocs)){
eventColors[traLocs]=rgb(t(col2rgb("green"))/255, alpha = alphas[traLocs])
eventScaledColors[traLocs]=rgb(t(col2rgb("green"))/255, alpha = scaledAlphas[traLocs])  
}
if(any(otherLocs)){
eventColors[otherLocs]=rgb(t(col2rgb("gray"))/255, alpha = alphas[otherLocs])
eventScaledColors[otherLocs]=rgb(t(col2rgb("gray"))/255, alpha = scaledAlphas[otherLocs])
}
if(any(telInsLocs)){
eventColors[telInsLocs]=rgb(t(col2rgb("purple"))/255, alpha = alphas[telInsLocs])
eventScaledColors[telInsLocs]=rgb(t(col2rgb("purple"))/255, alpha = scaledAlphas[telInsLocs])
}


pdf(scaledOutput)
circos.initializeWithIdeogram()
circos.par(points.overflow.warning=FALSE)
title(pid)
circos.genomicLink(coords1[largeEventSelection,], coords2[largeEventSelection,], col = eventScaledColors[largeEventSelection])
circos.genomicLink(coords1[mediumEventSelection,], coords2[mediumEventSelection,], col = eventScaledColors[mediumEventSelection],h=0.15)
circos.genomicLink(coords1[smallEventSelection,], coords2[smallEventSelection,], col = eventScaledColors[smallEventSelection],h=0.05)
dev.off()
circos.clear()  
