#!/usr/bin/env bash

set -e -o pipefail

BPS_OUT_TMP=${BPS_OUT}.tmp.gz

${SAMTOOLS_BINARY} view -F 0x600 -f 0x001 ${BAMFILE} \
    | ${MBUFFER_BINARY} -q -m 2G -l /dev/null \
    | ${SOPHIA_BINARY} \
        --medianisize ${medianIsize} \                          # dito (maybe metadata table, later)
        --stdisizepercentage ${stdIsizePercentage} \           # map to [tumor,control]StdISizePercentage
        --properpairratio ${properPairRatio} \                  # dito
        --defaultreadlength ${defaultReadLength} \
        --clipsize ${clipThreshold} \
        --basequality ${qualThreshold} \
        --basequalitylow ${qualThresholdLow} \
        --lowqualclipsize ${lowQualOverhangThreshold} \
        --isizesigma ${isizeSigmaThreshold} \
        --bpsupport ${bpSupportThreshold} \
    | ${GZIP_BINARY} --best > ${BPS_OUT}

mv ${BPS_OUT_TMP} ${BPS_OUT}
