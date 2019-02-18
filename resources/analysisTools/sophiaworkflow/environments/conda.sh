#!/usr/bin/env bash
# Copyright (c) 2019
# Distributed under the MIT License (license terms are at https://www.github.com/ODCF/SophiaWorkflow/LICENSE.txt).
# Load a Conda environment.

source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export SAMTOOLS_BINARY=samtools
export BEDTOOLS_BINARY=bedtools
export PYTHON_BINARY=python
export RSCRIPT_BINARY=Rscript
export SOPHIA_ANNOTATION_BINARY=sophiaAnnotate
export SOPHIA_BINARY=sophia
export GZIP_BINARY="gzip"
export MBUFFER_BINARY="mbuffer"
