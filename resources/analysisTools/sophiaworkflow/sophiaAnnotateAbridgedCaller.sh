#!/bin/bash

if [ "$LOAD_MODULE" == true ]
then
        module load $MODULE_ENV
        SAMTOOLS_BINARY=samtools
        SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
        RSCRIPT_BINARY=Rscript
        SAMTOOLS_BINARY=samtools
        BEDTOOLS_BINARY=bedtools
        PYTHON_BINARY=python3
fi

## ADDED TO TEST

source python3sourceme

# Get rid of the extension: .bed.gz This is done for most of the files in Roddy, but some are left here. We still
# need to figure out, if the files should all be passed as parameters.
tumorFileRaw=${tumorFile:0:${#tumorFile}-7}
# Get rid of the path without "filename"
#tumorFileRaw=${tumorFileRaw##*/}

#TEMPORARY FILES
ABRIDGED_ANNOTATION=${tumorFileRaw}_annotatedAbridged.bedpe
ABRIDGED_ANNOTATION_TEMP=${ABRIDGED_ANNOTATION}.tmp
ABRIDGED_ANNOTATION_CONTEXT=${tumorFileRaw}_annotatedAbridgedContext.bedpe
FILE_DUM=${tumorFileRaw}_temp
QC_JSON_FILE=${outputAnalysisBaseDirectory}/${sophiaOutputDirectory}/${JSON_PREFIX}qualitycontrol.json
QC_JSON_FILE_TMP=${QC_JSON_FILE}.tmp
#TEMPORARY FILES END

if [[ ! -e "${bloodFile}" ]]
then
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} \
                                    --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} \
                                    --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} \
                                    --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} \
                                    --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --germlinedblimit ${germlineDbLimit} > ${ABRIDGED_ANNOTATION_TEMP}
else
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} \
                                    --controlresults ${bloodFile} --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} \
	                            --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} \
	                            --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} \
	                            --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --defaultreadlengthcontrol ${controlDefaultReadLength} --germlinedblimit ${germlineDbLimit} > ${ABRIDGED_ANNOTATION_TEMP}
fi

set +e
grep $'^#' ${ABRIDGED_ANNOTATION_TEMP} > ${ABRIDGED_ANNOTATION}.WARNINGS
set -e

# CREATE THE START OF THE QC JSON FILE
touch $QC_JSON_FILE_TMP
echo -e "{\n  \"all\": {" > $QC_JSON_FILE_TMP
# FOREACH WARNING ADD A VALUE TO THE QC JSON FILE
controlMassiveInvPrefilteringLevel=0
tumorMassiveInvFilteringLevel=0
if [[ "`grep controlMassiveInvPrefilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | wc -l`" != "0" ]]
then
	controlMassiveInvPrefilteringLevel=`grep controlMassiveInvPrefilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | cut -f 2`
fi

if [[ "`grep tumorMassiveInvFilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | wc -l`" != "0" ]]
then
        tumorMassiveInvFilteringLevel=`grep tumorMassiveInvFilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | cut -f 2`
fi

echo "    \"controlMassiveInvPrefilteringLevel\": ${controlMassiveInvPrefilteringLevel}," >> $QC_JSON_FILE_TMP
echo "    \"tumorMassiveInvFilteringLevel\": ${tumorMassiveInvFilteringLevel}," >> $QC_JSON_FILE_TMP

if [[ "`cat ${ABRIDGED_ANNOTATION}.WARNINGS | wc -l`" == "0" ]]
then
	rm ${ABRIDGED_ANNOTATION}.WARNINGS
fi

grep -v $'^#' ${ABRIDGED_ANNOTATION_TEMP} | awk '($4 != "NA")' | grep -v "GL00" | sort -V -k 1,1 -k2,2 -k4,4 -k5,5 > ${ABRIDGED_ANNOTATION}
rm ${ABRIDGED_ANNOTATION_TEMP}

cut -f 1-3 ${ABRIDGED_ANNOTATION} | cat -n | sed 's/ //g' | awk -v OFS='\t'  '{print $2, $3, $4, $1}' > ${FILE_DUM}leftCoords
${BEDTOOLS_BINARY} intersect -a ${FILE_DUM}leftCoords -b ${intronExonRefBed} ${geneRefBedCancer} ${combinedSuperEnhancerRefBed} -loj  > ${FILE_DUM}directHits1Pre
${PYTHON_BINARY} ${TOOL_DIRECTHITCOLLAPSE_SCRIPT} ${FILE_DUM}directHits1Pre >${FILE_DUM}directHits1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCoding} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 > ${FILE_DUM}genesUpstream1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,10 > ${FILE_DUM}genesUpstreamCancer1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCoding} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 > ${FILE_DUM}genesDownstream1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,10 > ${FILE_DUM}genesDownstreamCancer1

cut -f 4-6 ${ABRIDGED_ANNOTATION} | cat -n | sed 's/ //g' | awk -v OFS='\t'  '{print $2, $3, $4, $1}' | sort -V -k1,1 -k2,2 -k3,3 -k4,4 > ${FILE_DUM}rightCoords
cut -f4 ${FILE_DUM}rightCoords > ${FILE_DUM}lineOrder
${BEDTOOLS_BINARY} intersect -a ${FILE_DUM}rightCoords -b ${intronExonRefBed} ${geneRefBedCancer} ${combinedSuperEnhancerRefBed} -loj > ${FILE_DUM}directHits2Pre
${PYTHON_BINARY} ${TOOL_DIRECTHITCOLLAPSE_SCRIPT} ${FILE_DUM}directHits2Pre >${FILE_DUM}directHits2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCoding} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 > ${FILE_DUM}genesUpstream2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,10 > ${FILE_DUM}genesUpstreamCancer2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCoding} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 > ${FILE_DUM}genesDownstream2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,10 > ${FILE_DUM}genesDownstreamCancer2

paste ${FILE_DUM}lineOrder <(cut -f1,2 ${FILE_DUM}directHits2) ${FILE_DUM}genesUpstream2 ${FILE_DUM}genesUpstreamCancer2 ${FILE_DUM}genesDownstream2 ${FILE_DUM}genesDownstreamCancer2 | sort -V -k1,1 | cut -f2- > ${FILE_DUM}right
cut -f3 ${FILE_DUM}directHits1 > ${FILE_DUM}leftSuperEnhancers
paste ${FILE_DUM}lineOrder <(cut -f3 ${FILE_DUM}directHits2) | sort -V -k1,1 | cut -f2- > ${FILE_DUM}rightSuperEnhancers
paste ${ABRIDGED_ANNOTATION} <(cut -f1,2 ${FILE_DUM}directHits1) ${FILE_DUM}genesUpstream1 ${FILE_DUM}genesUpstreamCancer1 ${FILE_DUM}genesDownstream1 ${FILE_DUM}genesDownstreamCancer1 ${FILE_DUM}right ${FILE_DUM}leftSuperEnhancers ${FILE_DUM}rightSuperEnhancers > ${FILE_DUM}tadPrep
${PYTHON_BINARY} ${TOOL_INTRACHROMOSOMALEVENTPICKER} ${FILE_DUM}tadPrep | ${BEDTOOLS_BINARY} intersect -a stdin -b ${roadmapEnhancerRefBed} -u | cut -f 4 > ${FILE_DUM}tadWhitelist
${PYTHON_BINARY} ${TOOL_TADANNOTATION_SCRIPT} ${FILE_DUM}tadPrep ${FILE_DUM}tadWhitelist ${consensusTadReferenceBed} ${smallEventThreshold} > ${FILE_DUM}tadAnnotations
paste ${FILE_DUM}tadPrep ${FILE_DUM}tadAnnotations > ${ABRIDGED_ANNOTATION_CONTEXT}

cat <(echo -e ${STANDARDHEADER})  <(cat ${ABRIDGED_ANNOTATION_CONTEXT}) | uniq | ${TOOL_FUSION_CANDIDATES_SCRIPT} > ${BEDPE_RESULT_FILE_FILTERED}
#cat <(echo -e ${STANDARDHEADER})  <(cat ${ABRIDGED_ANNOTATION_CONTEXT}) | uniq  > ${BEDPE_RESULT_FILE_FILTERED}
#DELETE TEMPORARY FILES
rm ${ABRIDGED_ANNOTATION}
rm ${ABRIDGED_ANNOTATION_CONTEXT}
rm ${FILE_DUM}*
#TEMPORARY FILES END

#TEMPORARY FILES QC
MERGED_DELEXTRACT=${BEDPE_RESULT_FILE_FILTERED}_delExtract
MERGED_DELEXTRACTINTERSECT=${BEDPE_RESULT_FILE_FILTERED}_delExtractIntersect
MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS=${BEDPE_RESULT_FILE_FILTERED}_rnaContaminationCandidates
#TEMPORARY FILES QC END
${PYTHON_BINARY} ${TOOL_RNACONT_DEL_EXTRACTOR_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED}  > ${MERGED_DELEXTRACT}
${BEDTOOLS_BINARY} intersect -a ${MERGED_DELEXTRACT} -b ${exonsOnlyPaddedRefBed} -wa -wb > ${MERGED_DELEXTRACTINTERSECT}
${PYTHON_BINARY} ${TOOL_RNACONT_DEL_COUNTER_SCRIPT} ${MERGED_DELEXTRACTINTERSECT} | grep -v $'\t1$' >  ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS}
MERGED_RNA_CONTAMINATED_GENES="`cat ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS} | wc -l`"


# ADD TO QC_JSON
rnaContaminatedGenesCount=0
rnaContaminatedGenesMoreThanTwoIntron=`sed 's|\#done||' ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS} | cut -f 1 | sed ':a;N;$!ba;s/\n/;/g'`
echo "    \"rnaContaminatedGenesMoreThanTwoIntron\": \"${rnaContaminatedGenesMoreThanTwoIntron}\"," >> $QC_JSON_FILE_TMP
MERGED_RNA_CONTAMINATED_GENES2=`echo $MERGED_RNA_CONTAMINATED_GENES | awk -F'\t' '{print \$1-1}' `
echo "    \"rnaContaminatedGenesCount\": ${MERGED_RNA_CONTAMINATED_GENES2}," >> $QC_JSON_FILE_TMP

if [[ "${MERGED_RNA_CONTAMINATED_GENES}" -ge "11" ]]
then
        echo "    \"rnaDecontaminationApplied\": true" >> $QC_JSON_FILE_TMP
	mv ${BEDPE_RESULT_FILE_FILTERED} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue
	${PYTHON_BINARY} ${TOOL_RNADECONT_STEP1_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue | ${BEDTOOLS_BINARY} intersect -f 0.9 -r -a stdin -b ${spliceJunctionsRefBed} -wa | cut -f4 | uniq > ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices
	${PYTHON_BINARY} ${TOOL_RNADECONT_STEP2_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices > ${BEDPE_RESULT_FILE_FILTERED}
	rm ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices
	echo "#${MERGED_RNA_CONTAMINATED_GENES} genes are affected by RNA contamination, RNA decontamination applied, expect losses and false positives!" >> ${ABRIDGED_ANNOTATION}.WARNINGS
else
        echo "    \"rnaDecontaminationApplied\": false" >> $QC_JSON_FILE_TMP
	rm ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS}
fi
#DELETE TEMPORARY FILES QC
rm ${MERGED_DELEXTRACT}
rm ${MERGED_DELEXTRACTINTERSECT}
#TEMPORARY FILES QC END

if [[ -e "${bloodFile}" ]]
then
	cat <(head -n 1 ${BEDPE_RESULT_FILE_FILTERED})  <(grep GERMLINE ${BEDPE_RESULT_FILE_FILTERED}) | uniq > ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}
	grep -v GERMLINE ${BEDPE_RESULT_FILE_FILTERED}  > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}

        # FILTERING FOR ACESEQ WILL COME HERE
        awk '$10>=3' ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ}
        # FILTERING FOR ACESEQ WILL END HERE

	set +e
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S (${analysisTag})" "3"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} "${pid} Somatic (${analysisTag})" "3"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_GERMLINE} "${pid} Germline (${analysisTag})" "3"
	pdfunite ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_scaled_merged.pdf
	rm ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_scaled.pdf
        pdfunite ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_unscaled_merged.pdf
	rm ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_unscaled.pdf
else
        # FILTERING FOR ACESEQ WILL COME HERE
        awk '$10>=3' ${BEDPE_RESULT_FILE_FILTERED} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ}
        # FILTERING FOR ACESEQ WILL END HERE

	set +e
	:
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S(${analysisTag})" "3"
        mv ${BEDPE_RESULT_FILE_FILTERED}_score_3_scaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_scaled_merged.pdf
	mv ${BEDPE_RESULT_FILE_FILTERED}_score_3_unscaled.pdf ${BEDPE_RESULT_FILE_FILTERED}_score_3_unscaled_merged.pdf
fi

# CLOSE QC_JSON_FILE
echo -e "  }\n}" >> $QC_JSON_FILE_TMP
mv  $QC_JSON_FILE_TMP  $QC_JSON_FILE
