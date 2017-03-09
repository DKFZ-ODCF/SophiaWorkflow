#!/usr/bin/env bash
source ${CONFIG_FILE}

${SAMTOOLS_BINARY} view -L ${mergedRef} -F 0x600 -f 0x001 ${BAMFILE} | \
    ${MBUFFER_BINARY} -q -m 1G -l /dev/null | \
    ${SOPHIA_BINARY} --mergedisizes ${mergedIsizes} --clipsize ${clipThreshold} \
                     --basequality ${qualThreshold} --basequalitylow ${qualThresholdLow} \
                     --lowqualclipsize ${lowQualOverhangThreshold} --isizesigma ${isizeSigmaThreshold} \
                     --bpsupport ${bpSupportThreshold} --defaultreadlength ${defaultReadLength} | \
    ${GZIP_BINARY} --best > ${FILENAME_OUT}
