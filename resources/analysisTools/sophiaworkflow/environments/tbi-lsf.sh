#!/usr/bin/env bash


module load "samtools/$SAMTOOLS_VERSION" || throw 100 "Could not load module '$MODULE_ENV'"
module load "SOPHIA/$SOPHIA_VERSION" || throw 100 "Could not load module '$MODULE_ENV'"
module load "bedtools/$BEDTOOLS_VERSION"|| throw 100 "Could not load module '$MODULE_ENV'"
module load "python/$PYTHON_VERSION" || throw 100 "Could not load module '$MODULE_ENV'"
module load "R/$R_VERSION" || throw 100 "Could not load module '$MODULE_ENV'"

export PYTHON_BINARY=python

export SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
export SOPHIA_BINARY=sophia

export SAMTOOLS_BINARY=samtools
export RSCRIPT_BINARY=Rscript
export SAMTOOLS_BINARY=samtools
export BEDTOOLS_BINARY=bedtools

export GZIP_BINARY="gzip"
export MBUFFER_BINARY="mbuffer"