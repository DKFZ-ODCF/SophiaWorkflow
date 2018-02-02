#!/usr/bin/env bash

module load $MODULE_ENV || throw 100 "Could not load module '$MODULE_ENV'"

export PYTHON_BINARY=python

export SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
export SOPHIA_BINARY=sophia

export SAMTOOLS_BINARY=samtools
export RSCRIPT_BINARY=Rscript
export SAMTOOLS_BINARY=samtools
export BEDTOOLS_BINARY=bedtools

export GZIP_BINARY="gzip"
export MBUFFER_BINARY="mbuffer"