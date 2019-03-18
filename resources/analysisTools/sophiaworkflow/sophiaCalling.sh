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
