#!/bin/bash

set -e -o pipefail

# Grep returns 1 if no match is found. Ensure that empty files are handled gracefully without exit due to `set -e`.
grepIgnoreEmpty() {
    (set +e    # The parentheses keep the `set -e` in the local subshell.
    grep "$@"
    if [[ $? -gt 1 ]]; then
        return $?
    else
        return 0
    fi)
}

#IMPORTANT FILES
BEDPE_RESULT_FILE_FILTERED="$BEDPE_RESULT_FILE_FILTERED.tmp"
BEDPE_RESULT_FILE_FILTERED_SOMATIC="$BEDPE_RESULT_FILE_FILTERED_SOMATIC.tmp"
BEDPE_RESULT_FILE_FILTERED_GERMLINE="$BEDPE_RESULT_FILE_FILTERED_GERMLINE.tmp"

BEDPE_RESULT_FILE_FILTERED_DEDUP="$BEDPE_RESULT_FILE_FILTERED_DEDUP.tmp"
BEDPE_RESULT_FILE_FILTERED_DEDUP_SOMATIC="$BEDPE_RESULT_FILE_FILTERED_DEDUP_SOMATIC.tmp"
BEDPE_RESULT_FILE_FILTERED_DEDUP_GERMLINE="$BEDPE_RESULT_FILE_FILTERED_DEDUP_GERMLINE.tmp"

BEDPE_RESULT_FILE_FILTERED_ACESEQ="$BEDPE_RESULT_FILE_FILTERED_ACESEQ.tmp"
BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ="$BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ.tmp"
BEDPE_RESULT_FILE_FILTERED_SOMATIC_OVERHANG_CANDIDATES="$BEDPE_RESULT_FILE_FILTERED_SOMATIC_OVERHANG_CANDIDATES.tmp"
#IMPORTANT FILES END

#TEMPORARY FILES
tumorFileRaw="$outputAnalysisBaseDirectory/$sophiaOutputDirectory/"$(basename "$tumorFile" .tsv.gz)
ABRIDGED_ANNOTATION=${tumorFileRaw}_annotatedAbridged.bedpe
ABRIDGED_ANNOTATION_CONTEXT="${tumorFileRaw}_annotatedAbridgedContext.bedpe.tmp"
FILE_DUM=${tumorFileRaw}_dum
#TEMPORARY FILES END

if [[ ! -e "${bloodFile}" ]]
then
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --germlinedblimit ${germlineDbLimit} > ${ABRIDGED_ANNOTATION}Pre
else
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} --controlresults ${bloodFile} --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --defaultreadlengthcontrol ${controlDefaultReadLength} --germlinedblimit ${germlineDbLimit} > ${ABRIDGED_ANNOTATION}Pre
fi


grepIgnoreEmpty $'^#' ${ABRIDGED_ANNOTATION}Pre | grepIgnoreEmpty -v $'^##' > ${ABRIDGED_ANNOTATION}.WARNINGS
grepIgnoreEmpty $'^##' ${ABRIDGED_ANNOTATION}Pre  | sed 's/^##//' > "$BEDPE_RESULT_FILE_FILTERED_SOMATIC_OVERHANG_CANDIDATES"


grepIgnoreEmpty -v $'^#' ${ABRIDGED_ANNOTATION}Pre | awk '($4 != "NA")' | awk '$10 > 0' | sort -V -k 1,1 -k2,2 -k4,4 -k5,5 > ${ABRIDGED_ANNOTATION}_preRemap


${PYTHON_BINARY} ${TOOL_DECOY_MAPPER_SCRIPT} ${decoyRangeRefBed} ${ABRIDGED_ANNOTATION}_preRemap | sort -V -k 1,1 -k2,2 -k4,4 -k5,5 > ${ABRIDGED_ANNOTATION}_tmp
rev ${ABRIDGED_ANNOTATION}_tmp | cut -f 1-6 | rev > ${ABRIDGED_ANNOTATION}_preRemapCoords
rev ${ABRIDGED_ANNOTATION}_tmp | cut -f 7- | rev > ${ABRIDGED_ANNOTATION}

rm ${ABRIDGED_ANNOTATION}Pre
rm ${ABRIDGED_ANNOTATION}_tmp
rm ${ABRIDGED_ANNOTATION}_preRemap


cut -f 1-3 ${ABRIDGED_ANNOTATION} | cat -n | sed 's/ //g' | awk -v OFS='\t'  '{print $2, $3, $4, $1}' > ${FILE_DUM}leftCoords
${BEDTOOLS_BINARY} intersect -a ${FILE_DUM}leftCoords -b ${intronExonRefBed} ${intronExonRefBedCancer} ${combinedSuperEnhancerRefBed} -loj  > ${FILE_DUM}directHits1Pre
${PYTHON_BINARY} ${TOOL_DIRECTHITCOLLAPSE_SCRIPT} ${FILE_DUM}directHits1Pre >${FILE_DUM}directHits1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesUpstream1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesUpstreamCancer1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesDownstream1
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}leftCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesDownstreamCancer1

cut -f 4-6 ${ABRIDGED_ANNOTATION} | cat -n | sed 's/ //g' | awk -v OFS='\t'  '{print $2, $3, $4, $1}' | sort -V -k1,1 -k2,2 -k3,3 -k4,4 > ${FILE_DUM}rightCoords
cut -f4 ${FILE_DUM}rightCoords > ${FILE_DUM}lineOrder
${BEDTOOLS_BINARY} intersect -a ${FILE_DUM}rightCoords -b ${intronExonRefBed} ${intronExonRefBedCancer} ${combinedSuperEnhancerRefBed} -loj > ${FILE_DUM}directHits2Pre
${PYTHON_BINARY} ${TOOL_DIRECTHITCOLLAPSE_SCRIPT} ${FILE_DUM}directHits2Pre >${FILE_DUM}directHits2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesUpstream2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesUpstreamCancer2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesDownstream2
${BEDTOOLS_BINARY} closest -a ${FILE_DUM}rightCoords -b ${geneRefBedCancer} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 8,9 | sed 's/^\t/.\t/' > ${FILE_DUM}genesDownstreamCancer2

paste ${FILE_DUM}lineOrder <(cut -f1,2 ${FILE_DUM}directHits2) ${FILE_DUM}genesUpstream2 ${FILE_DUM}genesUpstreamCancer2 ${FILE_DUM}genesDownstream2 ${FILE_DUM}genesDownstreamCancer2 | sort -V -k1,1 | cut -f2- > ${FILE_DUM}right
cut -f3 ${FILE_DUM}directHits1 > ${FILE_DUM}leftSuperEnhancers
paste ${FILE_DUM}lineOrder <(cut -f3 ${FILE_DUM}directHits2) | sort -V -k1,1 | cut -f2- > ${FILE_DUM}rightSuperEnhancers
paste ${ABRIDGED_ANNOTATION} <(cut -f1,2 ${FILE_DUM}directHits1) ${FILE_DUM}genesUpstream1 ${FILE_DUM}genesUpstreamCancer1 ${FILE_DUM}genesDownstream1 ${FILE_DUM}genesDownstreamCancer1 ${FILE_DUM}right ${FILE_DUM}leftSuperEnhancers ${FILE_DUM}rightSuperEnhancers > ${FILE_DUM}tadPrep

${PYTHON_BINARY} ${TOOL_INTRACHROMOSOMALEVENTPICKER} ${FILE_DUM}tadPrep | ${BEDTOOLS_BINARY} intersect -a stdin -b ${roadmapEnhancerRefBed} -u | cut -f 4 > ${FILE_DUM}tadWhitelist
${PYTHON_BINARY} ${TOOL_TADANNOTATION_SCRIPT} ${FILE_DUM}tadPrep ${FILE_DUM}tadWhitelist ${consensusTadReferenceBed} ${smallEventThreshold} > ${FILE_DUM}tadAnnotations
paste ${FILE_DUM}tadPrep ${FILE_DUM}tadAnnotations ${ABRIDGED_ANNOTATION}_preRemapCoords > ${ABRIDGED_ANNOTATION_CONTEXT}
rm ${ABRIDGED_ANNOTATION}_preRemapCoords

cat <(echo -e ${STANDARDHEADER})  <(cat ${ABRIDGED_ANNOTATION_CONTEXT}) | uniq | ${TOOL_FUSION_CANDIDATES_SCRIPT} > ${BEDPE_RESULT_FILE_FILTERED}

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
${PYTHON_BINARY} ${TOOL_RNACONT_DEL_COUNTER_SCRIPT} ${MERGED_DELEXTRACTINTERSECT} | grepIgnoreEmpty -v $'\t1$' >  ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS}
MERGED_RNA_CONTAMINATED_GENES="`grepIgnoreEmpty -v locus ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS} | wc -l`"
if [[ "${MERGED_RNA_CONTAMINATED_GENES}" -ge "11" ]]
then
	mv ${BEDPE_RESULT_FILE_FILTERED} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue
	${PYTHON_BINARY} ${TOOL_RNADECONT_STEP1_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue | ${BEDTOOLS_BINARY} intersect -f 0.9 -r -a stdin -b ${spliceJunctionsRefBed} -wa | cut -f4 | uniq > ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices
	${PYTHON_BINARY} ${TOOL_RNADECONT_STEP2_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED}.preRnaRescue ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices > ${BEDPE_RESULT_FILE_FILTERED}
	rm ${BEDPE_RESULT_FILE_FILTERED}.contaminatedIndices
	echo "#${MERGED_RNA_CONTAMINATED_GENES} are affected by RNA contamination, RNA decontamination applied, expect losses and false positives!" >> ${ABRIDGED_ANNOTATION}.WARNINGS
fi
#DELETE TEMPORARY FILES QC
rm ${MERGED_DELEXTRACT}
rm ${MERGED_DELEXTRACTINTERSECT}
#TEMPORARY FILES QC END

${PYTHON_BINARY} ${TOOL_DEDUP_RESULTS_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} ${tumorDefaultReadLength} > ${BEDPE_RESULT_FILE_FILTERED_DEDUP}

if [[ -e "${bloodFile}" ]]
then
	cat <(head -n 1 ${BEDPE_RESULT_FILE_FILTERED})  <(grepIgnoreEmpty GERMLINE ${BEDPE_RESULT_FILE_FILTERED}) | uniq > ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}
	grepIgnoreEmpty -v GERMLINE ${BEDPE_RESULT_FILE_FILTERED}  > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}
	
	cat <(head -n 1 ${BEDPE_RESULT_FILE_FILTERED_DEDUP})  <(grepIgnoreEmpty GERMLINE ${BEDPE_RESULT_FILE_FILTERED_DEDUP}) | uniq > ${BEDPE_RESULT_FILE_FILTERED_DEDUP_GERMLINE}
	grepIgnoreEmpty -v GERMLINE ${BEDPE_RESULT_FILE_FILTERED_DEDUP}  > ${BEDPE_RESULT_FILE_FILTERED_DEDUP_SOMATIC}
	
	awk '$10>2' ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ}
	set +e
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} "'${project}:${pid}(somatic)'" "$circlizeScoreThreshold"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_GERMLINE} "'${project}:${pid}(germline)'" "$circlizeScoreThreshold"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "'${project}:${pid}(g&s)'" "$circlizeScoreThreshold"
	set -e

	mv "$BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ" "${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ/.tmp/}"
    mv "$BEDPE_RESULT_FILE_FILTERED_SOMATIC" "${BEDPE_RESULT_FILE_FILTERED_SOMATIC/.tmp/}"
    mv "$BEDPE_RESULT_FILE_FILTERED_DEDUP_SOMATIC" "${BEDPE_RESULT_FILE_FILTERED_DEDUP_SOMATIC/.tmp/}"
    mv "$BEDPE_RESULT_FILE_FILTERED_GERMLINE" "${BEDPE_RESULT_FILE_FILTERED_GERMLINE/.tmp/}"
    mv "$BEDPE_RESULT_FILE_FILTERED_DEDUP_GERMLINE" "${BEDPE_RESULT_FILE_FILTERED_DEDUP_GERMLINE/.tmp/}"
    mv "$BEDPE_RESULT_FILE_FILTERED_SOMATIC_OVERHANG_CANDIDATES" "${BEDPE_RESULT_FILE_FILTERED_SOMATIC_OVERHANG_CANDIDATES/.tmp/}"

    pdfunite \
	    "${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_${circlizeScoreThreshold}_scaled.pdf" \
	    "${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_${circlizeScoreThreshold}_scaled.pdf" \
	    "${BEDPE_RESULT_FILE_FILTERED}_score_${circlizeScoreThreshold}_scaled.pdf" \
	    "${BEDPE_RESULT_FILE_FILTERED_PDF}"
	rm "${BEDPE_RESULT_FILE_FILTERED_SOMATIC}_score_${circlizeScoreThreshold}_scaled.pdf"
	rm "${BEDPE_RESULT_FILE_FILTERED_GERMLINE}_score_${circlizeScoreThreshold}_scaled.pdf"
    rm "${BEDPE_RESULT_FILE_FILTERED}_score_${circlizeScoreThreshold}_scaled.pdf"

else
	awk '$10>2' ${BEDPE_RESULT_FILE_FILTERED} > ${BEDPE_RESULT_FILE_FILTERED_ACESEQ}

	set +e
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "'${project}:${pid}(g&s)'" "$circlizeScoreThreshold"
	set -e

	mv "$BEDPE_RESULT_FILE_FILTERED_ACESEQ" "${BEDPE_RESULT_FILE_FILTERED_ACESEQ/.tmp/}"
    mv "${BEDPE_RESULT_FILE_FILTERED}_score_${circlizeScoreThreshold}_scaled.pdf" "$BEDPE_RESULT_FILE_FILTERED_PDF"
fi

mv "$BEDPE_RESULT_FILE_FILTERED" "${BEDPE_RESULT_FILE_FILTERED/.tmp/}"
mv "$BEDPE_RESULT_FILE_FILTERED_DEDUP" "${BEDPE_RESULT_FILE_FILTERED_DEDUP/.tmp/}"

############################################ CREATE THE QC JSON FILE
QC_JSON_FILE="$QC_JSON_FILE.tmp"
touch $QC_JSON_FILE
echo -e "{\n  \"all\": {" > $QC_JSON_FILE
# FOREACH WARNING ADD A VALUE TO THE QC JSON FILE
controlMassiveInvPrefilteringLevel=0
tumorMassiveInvFilteringLevel=0
if [[ "`grepIgnoreEmpty controlMassiveInvPrefilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | wc -l`" != "0" ]]
then
	controlMassiveInvPrefilteringLevel=`grepIgnoreEmpty controlMassiveInvPrefilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | cut -f 2`
fi

if [[ "`grepIgnoreEmpty tumorMassiveInvFilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | wc -l`" != "0" ]]
then
        tumorMassiveInvFilteringLevel=`grepIgnoreEmpty tumorMassiveInvFilteringLevel ${ABRIDGED_ANNOTATION}.WARNINGS | cut -f 2`
fi

echo "    \"controlMassiveInvPrefilteringLevel\": ${controlMassiveInvPrefilteringLevel}," >> $QC_JSON_FILE
echo "    \"tumorMassiveInvFilteringLevel\": ${tumorMassiveInvFilteringLevel}," >> $QC_JSON_FILE

# ADD TO QC_JSON
rnaContaminatedGenesMoreThanTwoIntron=`sed 's|\#done||' ${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS} | cut -f 1 | sed ':a;N;$!ba;s/\n/;/g'`
echo "    \"rnaContaminatedGenesMoreThanTwoIntron\": \"${rnaContaminatedGenesMoreThanTwoIntron}\"," >> $QC_JSON_FILE

MERGED_RNA_CONTAMINATED_GENES2=`echo $MERGED_RNA_CONTAMINATED_GENES | awk -F'\t' '{print \$1-1}' `
echo "    \"rnaContaminatedGenesCount\": ${MERGED_RNA_CONTAMINATED_GENES2}," >> $QC_JSON_FILE

if [[ "${MERGED_RNA_CONTAMINATED_GENES}" -ge "11" ]]
then
    echo "    \"rnaDecontaminationApplied\": true" >> $QC_JSON_FILE
else
    echo "    \"rnaDecontaminationApplied\": false" >> $QC_JSON_FILE
fi

# CLOSE QC_JSON_FILE
echo -e "  }\n}" >> $QC_JSON_FILE

mv "$QC_JSON_FILE" "${QC_JSON_FILE/.tmp/}"




# Remove temporary files
rm "${MERGED_RNACONTCANDIDATES_MORETHAN2INTRONS}"
