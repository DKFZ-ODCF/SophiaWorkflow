#!/usr/bin/env bash

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