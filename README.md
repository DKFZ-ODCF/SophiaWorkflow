# Description

[Sophia](https://bitbucket.org/utoprak/sophia) is a tool for identifying somatic structural variations in tumor-control sample pairs. It was developed by Umut Toprak at the Computational Oncology group, Theoretical Bioinformatics, German Cancer Research Center (DKFZ). This repository code allows to run Sophia with [Roddy](https://github.com/TheRoddyWMS/Roddy).

## Installation

### Software Stack & Reference Data

Please refer to the documentation of the [Sophia](https://bitbucket.org/utoprak/sophia) tool for information on the bioinformatic software stack and reference data.

Please refer to the documentation of [Roddy](https://github.com/TheRoddyWMS/Roddy) for software requirements of the workflow management system.

The workflow depends on the [COWorkflowBasePlugin](https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin). Please refer to the `buildinfo.txt` in the plugin root directory to see which version is required for the specific version of this plugin you want to install.

## Configuration

All output goes to `outputAnalysisBaseDirectory`, which is the output base directory (`outputBaseDirectory`) given on the command line via `--useiodir` with the dataset identifier as subdirectory, i.e.

```
outputAnalysisBaseDirectory=$outputBaseDirectory/$dataSet
```

The output includes the `roddyExecutionStore` with the execution metadata. Usually, one additionally configures a `sophiaOutputDirectory` into which all SOPHIA output datat is written. Thus the output will go to:

```
$outputBaseDirectory/$dataSet/roddyExecutionStore
$outputBaseDirectory/$dataSet/$sophiaOutputDirectory
```

If you want all output to go into the `outputAnalysisBaseDirectory` just set `sophiaOutputDirectory` to an empty string or ".".


### Run flags / switches / passable values

| Switch                     | Default Description |
|----------------------------|---------------------|
| REFERENCES_BASE            | Path to the installation of the reference files. |
| bamfile_list               | A semicolon separated list of bamfiles (1. control, 2. tumor, ...) |
| sample_list | Semicolon-separated list of sample names |
| possibleControlSampleNamePrefixes | Space-separated list of prefix identifying control samples. Used for matching sample names in files when retrieving BAM metadata from pathnames. Always required. |
| possibleTumorSampleNamePrefixes | Space-separated list of prefixes identifying tumor samples. Used for matching sample names in files when retrieving BAM metadata from pathnames. Always required. |
| controlDefaultReadLength    | Default read length \[base pairs\] |
| tumorDefaultReadLength      | Default read length \[base pairs\] |
| controlMedianIsize          | Median of control insert size distribution \[base pairs\] |
| tumorMedianIsize            | Median of tumor insert size distribution \[base pairs\] |
| controlStdIsizePercentage   | \[median insert size / standard deviation insert size * 100\] |   
| tumorStdIsizePercentage     | \[median insert size / standard deviation insert size * 100\] |
| controlProperPairPercentage | Properly paired read pairs in control \[#(properly paired) / #(paired) * 100\] |
| tumorProperPairPercentage   | Properly paired read pairs in tumor \[#(properly paired) / #(paired) * 100\]|
--------------------------------

## Starting the Workflow

Most metadata are best provided via dedicated variables. 

```bash
roddy.sh run $configName@$analysisName $pid \
  --useconfig=/path/to/your/applicationProperties.ini \
  --configurationValues="\
      REFERENCES_BASE:/the/path/to/your/sophia/35/reference/files,\
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

With version 2, it is not possible anymore to provide the insert sizes via the `insertsizesfile_list` variable, like it was for the version 1. We strongly suggest you configure the workflow metadata manually, like shown above. However, if you want to retrieve the BAM files and their metadata from the filesystem, you can set `extractSamplesFromOutputFiles` to "true". This mode is less safe and clear than the more explicit way of calling the workflow and it should only work smoothly with the output of the AlignmentAndQCWorkflows Roddy-plugin. Furthermore, the `alignmentOutputDirectory` (default "alignment"), `insertSizesOutputDirectory` (default "insertsize_distribution", and `qcOutputDirectory` (default "qualitycontrol") need to be set correctly (and be of type "string" as preconfigured in the XML).

## Change Log

* 2.2.1

  * Bugfix: intersectionCollapsing.py produces empty output instead of dot-line on empty input. Avoid downsteam int-cast error

* 2.2.0

  * Fix ignoring of `sophiaOutputDirectory` for `$sampleType_$pid_bps_annotatedAbridged.bedpe.WARNINGS` file resulting in wrong location of the file
  * Update to Roddy 3.5
  * Update to COWorkflowBasePlugin 1.4.0
  * Bugfix: grepIgnoreEmpty (set +e on global environment was removed)
  * Bugfix: correct location of tumor file in caller
  * Added executability check for sample metadata configuration values
  * Refactorings to improve Bash code robustness 
  * Consistent usage of grepIgnoreEmpty to handle input files without matches for low-coverage/WES input data
  * Update to newer workflow API (Roddy 3.5)
  * Documentation
  * GPL 2+

* Version 2.1.2

  * Bugfix: In no-control case a PDF file path was not defined.

* Version 2.1.1

  * Bugfix: Fixed another situation where PDFs were not generated

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
  
* Version 1.2.18

  * Consistent usage of grepIgnoreEmpty to handle input files without matches for low-coverage/WES input data
  
* Version 1.2.17

  * Fixed some greps to ignore non-matching input files for low-coverage/WES input data

* Version 1.2.16

  * Migration to Roddy 3 and its LSF support
  * Environment scripts for DKFZ-ODCF cluster
  
* Version 1.0.16

  * Sophia version 34.
  
# License

See the [licence file](LICENSE.md).
