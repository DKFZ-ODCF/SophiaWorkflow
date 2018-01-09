#!/usr/bin/env bash

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV

	export PYTHON_BINARY=python

	export SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
	export SOPHIA_BINARY=sophia

	export SAMTOOLS_BINARY=samtools
	export RSCRIPT_BINARY=Rscript
	export BEDTOOLS_BINARY=bedtools
else
    export PYTHON_BINARY="python3"

    export SOPHIA_BINARY="${PIPELINE_BASE:?No PIPELINE_BASE variable defined}/sophia"
    export SOPHIA_ANNOTATION_BINARY="${ANNOTATION_PIPELINE_BASE:?No ANNOTATION_PIPELINE_BASE variable defined}/sophiaAnnotate"

    export SAMTOOLS_BINARY="samtools-1.2"
    export RSCRIPT_BINARY="Rscript-3.0.0"
    export BEDTOOLS_BINARY="bedtools-2.24.0"
fi

export GZIP_BINARY="gzip"
export MBUFFER_BINARY="mbuffer"