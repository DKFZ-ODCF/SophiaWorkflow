package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.b080.co.common.WorkflowUsingMergedBams
import de.dkfz.b080.co.files.BasicBamFile
import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationValue
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
     * @param sampleName        tumor, control, ...
     * @return
     */
    static LinkedHashMap<String, String> getBamParameters(Configuration config, String sampleName) {
        // TODO: Make this workflow run with a metadata table. Probably the WorkflowUsingMergedBams has to be adapted (currently it only accepts a single pair!
        // if (Roddy.isMetadataCLOptionSet())
        //    SophiaMetadataTable metadataTable = new SophiaMetadataTable(MetadataTableFactory.getTable(context.analysis))

        LinkedHashMap<String, String> result = new LinkedHashMap()
        ["defaultReadLength",
         "medianIsize",
         "stdIsizePercentage",
         "properPairRatio"
        ].collectEntries(result) { varname ->
            new MapEntry(varname, config.configurationValues.get("${sampleName}${varname.capitalize()}").toString())
        }
        return result
    }

    @Override
    boolean execute(ExecutionContext context, BasicBamFile bamControlMerged, BasicBamFile bamTumorMerged) {
        boolean fullWorkflow = bamControlMerged

        context.configuration.configurationValues << new ConfigurationValue("tumorSample", (bamTumorMerged.fileStage as COFileStageSettings).sample.name)

        if (fullWorkflow) { // Control and tumor
            context.configuration.configurationValues << new ConfigurationValue(ANALYSIS_TAG, "tumorControl")
            context.configuration.configurationValues << new ConfigurationValue("controlSample", (bamControlMerged.fileStage as COFileStageSettings).sample.name)

            BaseFile sophiaControlFile =
                    run(TOOL_SOPHIA, bamControlMerged, getBamParameters(context.configuration, "control")) as BaseFile
            BaseFile sophiaTumorFile =
                    run(TOOL_SOPHIA, bamTumorMerged, getBamParameters(context.configuration, "tumor")) as BaseFile

            run("sophiaAnnotator", sophiaControlFile, sophiaTumorFile)
        } else { //NoControl!
            context.getConfiguration().getConfigurationValues().add(new ConfigurationValue(ANALYSIS_TAG, "tumorOnly"))
            BaseFile sophiaTumorFile =
                    run(TOOL_SOPHIA, bamTumorMerged, getBamParameters(context.configuration, "tumor")) as BaseFile

            run("sophiaAnnotatorNoControl", sophiaTumorFile)
        }
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
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("${BAMFILE_LIST} is not set."))
            return false
        }

        return true
    }
}
