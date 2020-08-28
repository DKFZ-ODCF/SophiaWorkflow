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

mv "$BPS_OUT_TMP" "$BPS_OUT" || throw 100 "Could not move file: '$BPS_OUT_TMP'"

outputSize=$(gunzip -c "$BPS_OUT" | head -n 10 | wc -l)
if [[ $outputSize -eq 0 ]]; then
    # The file should have a header.
    echo "sophia binary may have failed. No output whatsoever" >> /dev/stderr
    exit 10
elif [[ $outputSize -eq 1 ]]; then
    # The file must have a header.
    # Empty output means, the data may have been of too low quality or depth.
    # Instead of continuing and failing in a later job, rather fail here early with a dedicated exit code.
    # BUT, OTP is not prepared to handle this condition, so for the time being return exit 0 and continue!
    echo "sophia binary call yielded empty results set" >> /dev/stderr
    exit 0
fi

