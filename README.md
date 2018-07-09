# Description

Sophia is a tool for identifying structural variations. It was developed by Umut Toprak at the Computational Oncology group, Theoretical Bioinformatics, German Cancer Research Center (DKFZ).

## Run flags / switches / passable values

| Switch                     | Default Description |
|----------------------------|---------------------|
| bamfile_list               | A semicolon separated list of bamfiles (1. control, 2. tumor, ...) |
| possibleControlSampleNamePrefixes | Space-separated list of prefix identifying control samples |
| possibleTumorSampleNamePrefixes | Space-separated list of prefix identifying tumor samples |
| sample_list | Semicolon-separated list of sample names |
| controlDefaultReadLength   | Default read length |
| tumorDefaultReadLength     | Default read length |
| controlMedianIsize         | Median of control insert size distribution |
| tumorMedianIsize           | Median of tumor insert size distribution |
| controlStdIsizePercentage  | Quotient of median and standard deviation of the control insert size distribution |   
| tumorStdIsizePercentage    | Quotient of median and standard deviation of the tumor insert size distribution |
| controlProperPairRatio     | Proportion of properly paired read pairs in control |
| tumorProperPairRatio       | Proportion of properly paired read pairs in tumor |
--------------------------------

## Starting the Workflow

Most metadata are best provided via dedicated variables. 

```bash
roddy.sh run $configName@$analysisName $pid \
  --useconfig=/path/to/your/applicationProperties.ini \
  --configurationValues="\
      bamfile_list:$controlBam;$tumorBam,\
      controlDefaultReadLength:${controlDefaultReadLength},\
      tumorDefaultReadLength:${tumorDefaultReadLength},\
      controlMedianIsize:${controlMedianIsize},\
      tumorMedianIsize:${tumorMedianIsize},\
      controlStdIsizePercentage:${controlStdIsizePercentage},\
      tumorStdIsizePercentage:${tumorStdIsizePercentage},\
      controlProperPairRatio:${controlProperPairRatio},\
      tumorProperPairRatio:${tumorProperPairRatio}"

```

It is not possible anymore to provide the insert sizes via the `insertsizesfile_list` variable, like it was for the earlier versions.

If you want to retrieve the BAM files and their metadata from the filesystem, you can also use the following variables:

  * `extractSamplesFromOutputFiles` needs to be set to "true"
  * `possibleControlSampleNamePrefixes`
  * `possibleTumorSampleNamePrefixes`
  * `sample_list`

However, this mode is less safe and clear than the previous more explicit way of calling the workflow. On the longer run, we will implemented the
metadata provisioning via metadata table which is more convenient for calling and less error prone.

## Changelist

* Version 2.0.0

  * Sophia update to 35
  
* Version 1.2.