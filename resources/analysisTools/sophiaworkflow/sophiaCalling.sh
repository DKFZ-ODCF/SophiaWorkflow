#!/usr/bin/env bash
#
# Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

set -e -o pipefail -u

BPS_OUT_TMP="$BPS_OUT.tmp.gz"

ignoreEc() {
  local observedCode="$?"
  local ignorableErrorCode="${1:?No error code given}"
  if [[ $observedCode -eq $ignorableErrorCode ]]; then
    return 0
  else
    return $observedCode
  fi
}

"$SAMTOOLS_BINARY" view -F 0x600 -f 0x001 "$BAMFILE" \
    | "$MBUFFER_BINARY" -q -m 2G -l /dev/null \
    | "$SOPHIA_BINARY" \
        --medianisize "$medianIsize" \
        --stdisizepercentage "$stdIsizePercentage" \
        --properpairpercentage "$properPairPercentage" \
        --defaultreadlength "$defaultReadLength" \
        --clipsize "$clipThreshold" \
        --basequality "$qualThreshold" \
        --basequalitylow "$qualThresholdLow" \
        --lowqualclipsize "$lowQualOverhangThreshold" \
        --isizesigma "$isizeSigmaThreshold" \
        --bpsupport "$bpSupportThreshold" \
    | "$GZIP_BINARY" --best > "$BPS_OUT_TMP"

# gunzip or zcat exit with code (141), because `head` closes the pipe.
BROKEN_PIPE=141
outputSize=$(gunzip -c "$BPS_OUT_TMP" | head -n 10 | wc -l || ignoreEc $BROKEN_PIPE)
if [[ "$outputSize" -eq 0 ]]; then
    # The file should have at least a header. Otherwise this is an indication that ...
    echo "Sophia binary may have failed. No output whatsoever" >> /dev/stderr
    exit 10
elif [[ "$outputSize" -eq 1 ]]; then
    # If the file has a header but otherwise no output, the data may have been of too low quality or depth.
    # Instead of continuing and failing in a later job we could fail here early with a dedicated exit code.
    # BUT, OTP is not prepared to handle this condition, so for the time being, return exit 0 and continue!
    echo "Sophia binary call yielded empty results set" >> /dev/stderr
    mv "$BPS_OUT_TMP" "$BPS_OUT" || throw 100 "Could not move file: '$BPS_OUT_TMP'"
    exit 0
else
    mv "$BPS_OUT_TMP" "$BPS_OUT" || throw 100 "Could not move file: '$BPS_OUT_TMP'"
fi

