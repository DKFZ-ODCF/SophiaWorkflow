#!/bin/bash
source ${CONFIG_FILE}

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

source python3sourceme

# Get rid of the extension: .bed.gz This is done for most of the files in Roddy, but some are left here. We still
# need to figure out, if the files should all be passed as parameters.
tumorFileRaw=${tumorFile:0:${#tumorFile}-7}
# Get rid of the path without "filename"
tumorFileRaw=${tumorFileRaw##*/}

#TEMPORARY FILES
ABRIDGED_ANNOTATION=${tumorFileRaw}_annotatedAbridged.bedpe
ABRIDGED_ANNOTATION_CONTEXT=${tumorFileRaw}_annotatedAbridgedContext.bedpe
FILE_LEFTCONTEXT=${tumorFileRaw}_leftContext
FILE_LEFTCONTEXT_PRE=${tumorFileRaw}_leftContextPre
FILE_RIGHTCONTEXT=${tumorFileRaw}_rightContext
FILE_RIGHTCONTEXT_PRE=${tumorFileRaw}_rightContextPre
FILE_DUMLEFT=${tumorFileRaw}_dumLeft
FILE_DUMRIGHT=${tumorFileRaw}_dumRight
#TEMPORARY FILES END

if [[ ! -e "${bloodFile}" ]]
then
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} \
	                            --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} \
	                            --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} \
	                            --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} \
	                            --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --germlinedblimit ${germlineDbLimit} \
	                              | awk '($4 != "NA")' | grep -v "GL00" | sort -V -k 1,1 -k2,2 -k4,4 -k5,5 > ${ABRIDGED_ANNOTATION}
else
	${SOPHIA_ANNOTATION_BINARY} --tumorresults ${tumorFile} --mref ${mRef} \
	                            --controlresults ${bloodFile} --pidsinmref ${pidsInMref} --bpfreq ${bpFreq} \
	                            --artifactlofreq ${artifactLoFreq} --artifacthifreq ${artifactHiFreq} \
	                            --clonalitystrictlofreq ${clonalityStrictLoFreq} --clonalitylofreq ${clonalityLoFreq} --clonalityhifreq ${clonalityHiFreq} \
	                            --germlineoffset ${germlineFuzziness} --defaultreadlengthtumor ${tumorDefaultReadLength} --defaultreadlengthcontrol ${controlDefaultReadLength} --germlinedblimit ${germlineDbLimit} \
	                              | awk '($4 != "NA")' | grep -v "GL00" | sort -V -k 1,1 -k2,2 -k4,4 -k5,5 > ${ABRIDGED_ANNOTATION}
fi

cut -f 1-3 ${ABRIDGED_ANNOTATION} | cat -n | sed 's/ //g' \
                | awk -v OFS='\t'  '{print $2, $3, $4, $1}' \
                | ${BEDTOOLS_BINARY} intersect -a stdin -b ${intronExonRefBed} -sorted -loj -g ${chromSizesRef} | tr -s '\t' '\t' | cut -f 1-4,8 > ${FILE_DUMLEFT}
                
${BEDTOOLS_BINARY} merge -d -1 -i ${FILE_DUMLEFT} -c 4,5 -o collapse -delim "," \
                | ${BEDTOOLS_BINARY} closest -a stdin -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 1-5,9-10 \
                | ${BEDTOOLS_BINARY} closest -a stdin -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 1-7,11-12 > ${FILE_LEFTCONTEXT_PRE}
                
${PYTHON_BINARY} ${TOOL_SOPHIA_ANNOTATE_COORDINATE_CORRECTION} ${FILE_LEFTCONTEXT_PRE} > ${FILE_LEFTCONTEXT}

cut -f 4-6 ${ABRIDGED_ANNOTATION}  | cat -n | sed 's/ //g' | awk -v OFS='\t'  '{print $2, $3, $4, $1}' | sort -V -k1,1 -k2,2 -k3,3 \
                | ${BEDTOOLS_BINARY} intersect -a stdin -b ${intronExonRefBed} -sorted -loj -g ${chromSizesRef} | tr -s '\t' '\t'  | cut -f 1-4,8 > ${FILE_DUMRIGHT}
                
${BEDTOOLS_BINARY} merge -d -1 -i ${FILE_DUMRIGHT} -c 4,5 -o collapse -delim "," \
                | ${BEDTOOLS_BINARY} closest -a stdin -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -id -t last -k 1 | cut -f 1-5,9-10 \
                | ${BEDTOOLS_BINARY} closest -a stdin -b ${geneRefBed} -g ${chromSizesRef} -io -D ref -iu -t first -k 1 | cut -f 1-7,11-12 > ${FILE_RIGHTCONTEXT_PRE}
                
${PYTHON_BINARY} ${TOOL_SOPHIA_ANNOTATE_COORDINATE_CORRECTION} ${FILE_RIGHTCONTEXT_PRE} | sort -n -k 4,4 > ${FILE_RIGHTCONTEXT}

paste ${ABRIDGED_ANNOTATION} <(cut -f 5- ${FILE_LEFTCONTEXT}) <(cut -f 5- ${FILE_RIGHTCONTEXT}) > ${ABRIDGED_ANNOTATION_CONTEXT}
cat <(echo -e ${STANDARDHEADER})  <(cat ${ABRIDGED_ANNOTATION_CONTEXT}) | uniq > ${BEDPE_RESULT_FILE_FILTERED}

#DELETE TEMPORARY FILES
rm ${ABRIDGED_ANNOTATION}
rm ${ABRIDGED_ANNOTATION_CONTEXT}
rm ${FILE_LEFTCONTEXT}
rm ${FILE_LEFTCONTEXT_PRE}
rm ${FILE_RIGHTCONTEXT}
rm ${FILE_RIGHTCONTEXT_PRE}
rm ${FILE_DUMLEFT}
rm ${FILE_DUMRIGHT}
#TEMPORARY FILES END


if [[ -e "${bloodFile}" ]]
then
	cat <(echo -e ${STANDARDHEADER})  <(grep GERMLINE ${BEDPE_RESULT_FILE_FILTERED}) | uniq > ${BEDPE_RESULT_FILE_FILTERED_GERMLINE}
	grep -v GERMLINE ${BEDPE_RESULT_FILE_FILTERED} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC}
        # FILTERING FOR ACESEQ WILL COME HERE
        awk '$10>=3' ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ}
        # FILTERING FOR ACESEQ WILL END HERE
	set +e
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S(${analysisTag}) ES:1" "1"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S(${analysisTag}) ES:3" "3"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} "${pid} Somatic(${analysisTag}) ES:1" "1"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} "${pid} Somatic(${analysisTag}) ES:3" "3"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_GERMLINE} "${pid} Germline(${analysisTag}) ES:1" "1"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED_GERMLINE} "${pid} Germline(${analysisTag}) ES:3" "3"
	# QC WILL COME HERE
        touch $QC_JSON_FILE
	# QC WILL END HERE
else
        # FILTERING FOR ACESEQ WILL COME HERE
        awk '$10>=3' ${BEDPE_RESULT_FILE_FILTERED_SOMATIC} > ${BEDPE_RESULT_FILE_FILTERED_SOMATIC_ACESEQ}
        # FILTERING FOR ACESEQ WILL END HERE
	set +e
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S(${analysisTag}) ES:1" "1"
	${RSCRIPT_BINARY} ${TOOL_CIRCLIZE_SCRIPT} ${BEDPE_RESULT_FILE_FILTERED} "${pid} G&S(${analysisTag}) ES:3" "3"
        # NO CONTROL QC WILL COME HERE
        touch QC_JSON_FILE
        # NO CONTORL QC WILL END HERE
fi
