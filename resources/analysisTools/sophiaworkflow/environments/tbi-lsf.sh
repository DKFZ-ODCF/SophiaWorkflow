#!/usr/bin/env bash
#
# Copyright (C) 2018 Michael Heinold, Philip R. Kensche and DKFZ Heidelberg
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

module load "samtools/$SAMTOOLS_VERSION" || throw 100 "Could not load module 'samtools'"
export SAMTOOLS_BINARY=samtools

module load "bedtools/$BEDTOOLS_VERSION"|| throw 100 "Could not load module 'bedtools'"
export BEDTOOLS_BINARY=bedtools

module load "python/$PYTHON_VERSION" || throw 100 "Could not load module 'python'"
export PYTHON_BINARY=python

module load "R/$R_VERSION" || throw 100 "Could not load module 'R'"
export RSCRIPT_BINARY=Rscript

module load "sophia/$SOPHIA_VERSION" || throw 100 "Could not load module 'SOPHIA'"
export SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
export SOPHIA_BINARY=sophia

export GZIP_BINARY="gzip"
export MBUFFER_BINARY="mbuffer"