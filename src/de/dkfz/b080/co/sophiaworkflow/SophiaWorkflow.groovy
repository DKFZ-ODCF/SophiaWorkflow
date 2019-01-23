package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.b080.co.common.WorkflowUsingMergedBams
import de.dkfz.b080.co.files.BasicBamFile
import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationError
import de.dkfz.roddy.config.ConfigurationValue
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.knowledge.files.BaseFile
import groovy.transform.CompileStatic

/**
 * Sophia workflow
 */
@CompileStatic
class SophiaWorkflow extends WorkflowUsingMergedBams {

    public static final String TOOL_SOPHIA = "sophia"
    public static final String ANALYSIS_TAG = "analysisTag"
    public static final String BAMFILE_LIST = "bamfile_list"

    /**
     * Compose a LinkedHashMap with parameter values for (e.g.) tumor and control taken from the configuration.
     * This simplifies the collection of configuration values with "tumor" and "control" prefixes. On the long
     * run this should be substituted by getting the data from a metadata table. Default values are not implemented.
     *
     * @param context
     * @param sampleType        tumor, control, ...
     * @return
     */
    static LinkedHashMap<String, String> getBamParameters(Configuration config, Sample.SampleType sampleType) {
        LinkedHashMap<String, String> result = new LinkedHashMap()
        ['defaultReadLength',
         'medianIsize',
         'stdIsizePercentage',
         'properPairPercentage'
        ].collectEntries(result) { varname ->
            String cvalueName = "${sampleType.name().toLowerCase()}${varname.capitalize()}"
            String cvalue = config.configurationValues.getString(cvalueName)
            if (cvalue.empty)
                throw new ConfigurationError("Configuration value '${cvalueName}' is unset", config)
            new MapEntry(varname, cvalue)
        }
        result
    }

    @Override
    boolean execute(ExecutionContext context, BasicBamFile bamControlMerged, BasicBamFile bamTumorMerged) {
        boolean fullWorkflow = bamControlMerged

        Configuration configuration = context.configuration
        RecursiveOverridableMapContainerForConfigurationValues configurationValues = configuration.configurationValues

        configurationValues <<
                new ConfigurationValue('tumorSample', (bamTumorMerged.fileStage as COFileStageSettings).sample.name)

        if (fullWorkflow) { // Control and tumor
            configurationValues <<
                    new ConfigurationValue(ANALYSIS_TAG, 'tumorControl')
            configurationValues <<
                    new ConfigurationValue('controlSample', (bamControlMerged.fileStage as COFileStageSettings).sample.name)

            BaseFile sophiaControlFile =
                    run(TOOL_SOPHIA, bamControlMerged,
                            getBamParameters(configuration, bamControlMerged.sample.sampleType)) as BaseFile
            BaseFile sophiaTumorFile =
                    run(TOOL_SOPHIA, bamTumorMerged,
                            getBamParameters(configuration, bamTumorMerged.sample.sampleType)) as BaseFile

            run('sophiaAnnotator', sophiaControlFile, sophiaTumorFile)
        } else { // NoControl!
            configurationValues <<
                    new ConfigurationValue(ANALYSIS_TAG, 'tumorOnly')
            BaseFile sophiaTumorFile =
                    run(TOOL_SOPHIA, bamTumorMerged,
                            getBamParameters(configuration, bamTumorMerged.sample.sampleType)) as BaseFile

            run('sophiaAnnotatorNoControl', sophiaTumorFile)
        }
    }

    /* Move the following three functions to a superclass. */
    List<BasicBamFile> getControlBamFiles() {
        loadInitialBamFilesForDataset(context).findAll { it.sample.sampleType == Sample.SampleType.CONTROL } as List
    }

    List<BasicBamFile> getTumorBamFiles() {
        loadInitialBamFilesForDataset(context).findAll { it.sample.sampleType == Sample.SampleType.TUMOR } as List
    }

    List<BasicBamFile> getUnkownBamFiles() {
        loadInitialBamFilesForDataset(context).findAll { it.sample.sampleType == Sample.SampleType.UNKNOWN } as List
    }


    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean initialBamCheck = super.checkExecutability(context)
        if (!initialBamCheck)
            return false

        loadInitialBamFilesForDataset(context)

        def configurationValues = context.getConfigurationValues()
        boolean bamListIsSet = configurationValues.hasValue(BAMFILE_LIST)
        if(! bamListIsSet) {
            context.addError(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("${BAMFILE_LIST} is not set."))
            return false
        }

        try {
            getBamParameters(context.configuration, Sample.SampleType.TUMOR)
        } catch (ConfigurationError e) {
            context.addError(ExecutionContextError.EXECUTION_SETUP_INVALID.expand(e.message))
        }
        if (!controlBamFiles.empty) {
            try {
                getBamParameters(context.configuration, Sample.SampleType.CONTROL)
            } catch (ConfigurationError e) {
                context.addError(ExecutionContextError.EXECUTION_SETUP_INVALID.expand(e.message))
            }
        }

        return true
    }
}
