#!/usr/bin/env bash

#IMPORTANT FILES
#${BPS_OUT}
#IMPORTANT FILES END

source ${CONFIG_FILE}

set -e

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
					 | ${GZIP_BINARY} --best > ${BPS_OUT}
