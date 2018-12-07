# Description

[Sophia](https://bitbucket.org/utoprak/sophia) is a tool for identifying somatic structural variations in tumor-control sample pairs. It was developed by Umut Toprak at the Computational Oncology group, Theoretical Bioinformatics, German Cancer Research Center (DKFZ). This repository code allows to run Sophia with [Roddy](https://github.com/TheRoddyWMS/Roddy).

## Installation

### Software Stack

TBD

### Reference Data

TBD

## Configuration

### Run flags / switches / passable values

| Switch                     | Default Description |
|----------------------------|---------------------|
| bamfile_list               | A semicolon separated list of bamfiles (1. control, 2. tumor, ...) |
| sample_list | Semicolon-separated list of sample names |
| possibleControlSampleNamePrefixes | Space-separated list of prefix identifying control samples. Used for matching sample names in files when retrieving BAM metadata from pathnames. Always required. |
| possibleTumorSampleNamePrefixes | Space-separated list of prefixes identifying tumor samples. Used for matching sample names in files when retrieving BAM metadata from pathnames. Always required. |
| controlDefaultReadLength    | Default read length |
| tumorDefaultReadLength      | Default read length |
| controlMedianIsize          | Median of control insert size distribution |
| tumorMedianIsize            | Median of tumor insert size distribution |
| controlStdIsizePercentage   | Quotient of median and standard deviation of the control insert size distribution |   
| tumorStdIsizePercentage     | Quotient of median and standard deviation of the tumor insert size distribution |
| controlProperPairPercentage | Percentage of properly paired read pairs in control |
| tumorProperPairPercentage   | Percentage of properly paired read pairs in tumor |
--------------------------------

## Starting the Workflow

Most metadata are best provided via dedicated variables. 

```bash
roddy.sh run $configName@$analysisName $pid \
  --useconfig=/path/to/your/applicationProperties.ini \
  --configurationValues="\
      bamfile_list:$controlBam;$tumorBam,\
      sample_list:$controlName;$tumorName,
      possibleControlSampleNamePrefixes:$controlName,\
      possibleTumorSampleNamePrefixes:$tumorName,\
      controlDefaultReadLength:${controlDefaultReadLength},\
      tumorDefaultReadLength:${tumorDefaultReadLength},\
      controlMedianIsize:${controlMedianIsize},\
      tumorMedianIsize:${tumorMedianIsize},\
      controlStdIsizePercentage:${controlStdIsizePercentage},\
      tumorStdIsizePercentage:${tumorStdIsizePercentage},\
      controlProperPairPercentage:${controlProperPairPercentage},\
      tumorProperPairPercentage:${tumorProperPairPercentage}"

```

It is not possible anymore to provide the insert sizes via the `insertsizesfile_list` variable, like it was for the version 1.

If you want to retrieve the BAM files and their metadata from the filesystem, you can also set `extractSamplesFromOutputFiles` to "true". Note however that this mode is less safe and clear than the more explicit way of calling the workflow. 

## Changelist

* Version 2.1.0

  * Updated to COWorkflowBasePlugin 1.3.0 that provides an alternative sample-name extraction algorithm.

* Version 2.0.3

  * Bugfix: PDF generation in very low data WES failed. Added "no data" PDFs. 

* Version 2.0.2

  * Bugfix: PDF merging for samples without called germline variations.
  * Bugfix: Undefined `newSize` in decoyMapper.py on some datasets. 

* Version 2.0.1

  * Produce again a single merged PDF instead of three independent PDFs.

* Version 2.0.0

  * Sophia version 35.
  
* Version 1.2

  * Sophia version 34.