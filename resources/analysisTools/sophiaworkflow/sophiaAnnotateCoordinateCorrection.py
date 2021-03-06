#!/usr/bin/env python
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

from sys import argv

rightSideFile = argv[1]

with open(rightSideFile) as f:
    for line in f:
        lineChunks = line.rstrip().split('\t')
        geneMatches = ",".join(sorted(list(set(lineChunks[4].split(',')))))
        lineIndices = set(lineChunks[3].split(','))
        for index in sorted(list(lineIndices)):
            print("\t".join(lineChunks[0:3] + [index, geneMatches] + lineChunks[5:]))
