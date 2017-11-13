#!/usr/bin/env bash

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

set -e

BPS_OUT_TMP=${BPS_OUT}.tmp.gz

${SAMTOOLS_BINARY} view -L ${mergedRef} -F 0x600 -f 0x001 ${BAMFILE} \
  | ${MBUFFER_BINARY} -q -m 1G -l /dev/null \
  | ${SOPHIA_BINARY} --mergedisizes ${mergedIsizes} \
           --clipsize ${clipThreshold} \
           --basequality ${qualThreshold} \
           --basequalitylow ${qualThresholdLow} \
           --lowqualclipsize ${lowQualOverhangThreshold} \
           --isizesigma ${isizeSigmaThreshold} \
           --bpsupport ${bpSupportThreshold} \
           --defaultreadlength ${defaultReadLength} \
					 | ${GZIP_BINARY} --best > ${BPS_OUT_TMP}

mv ${BPS_OUT_TMP} ${BPS_OUT}
