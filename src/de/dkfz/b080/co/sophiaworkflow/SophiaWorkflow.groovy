package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.b080.co.common.WorkflowUsingMergedBams
import de.dkfz.b080.co.files.BasicBamFile
import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.config.ConfigurationValue
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.files.FileObject
import de.dkfz.roddy.plugins.LibrariesFactory
import groovy.transform.CompileStatic

/**
 * Sophia workflow
 */
@CompileStatic
class SophiaWorkflow extends WorkflowUsingMergedBams {

    public static final String TOOL_SOPHIA = "sophia"
    public static final String ANALYSIS_TAG = "analysisTag"
    public static final String ISIZES_LIST = "insertsizesfile_list"
    public static final String BAMFILE_LIST = "bamfile_list"

    /**
     * Load (create) the class for the insert sizes file.
     * This class actually exists in AlignmentAndQCWorkflows, but we do not want to add unnecessary dependencies.
     * Just use Roddys synthetic class capabilities here.
     * @return
     */
    Class<BaseFile> loadISizesFileClass() {
        LibrariesFactory.getInstance().loadRealOrSyntheticClass("InsertSizesValueFile", BaseFile.class as Class<FileObject>)
    }

    /**
     * Load the insert sizes file for a given bam file.
     * indexInBamList is only used, if the isizes files are passed as a cvalue
     *
     * @return
     */
    BaseFile getISizesForBam(ExecutionContext context, BasicBamFile bam, int indexInBamList) {
        def configurationValues = context.getConfigurationValues()
        boolean isizesListIsSet = configurationValues.hasValue(ISIZES_LIST);
        if (isizesListIsSet) {
            List<String> isizesList = configurationValues.getString(ISIZES_LIST).split(StringConstants.SPLIT_SEMICOLON) as List<String>;
            return BaseFile.constructSourceFile(loadISizesFileClass(), new File(isizesList[indexInBamList]), context, bam.fileStage)
        } else {
            return BaseFile.constructManual(loadISizesFileClass(), new BasicBamFile(bam)) as BaseFile
        }
    }

    @Override
    boolean execute(ExecutionContext context, BasicBamFile bamControlMerged, BasicBamFile bamTumorMerged) {
        boolean fullWorkflow = bamControlMerged

        context.configuration.configurationValues << new ConfigurationValue("tumorSample", (bamTumorMerged.fileStage as COFileStageSettings).sample.name)

        String tumorDefaultReadLength = context.configuration.configurationValues.get("tumorDefaultReadLength")
        String controlDefaultReadLength= context.configuration.configurationValues.get("controlDefaultReadLength")


        if (fullWorkflow) { // Control and tumor
            context.configuration.configurationValues << new ConfigurationValue(ANALYSIS_TAG, "tumorControl")
            context.configuration.configurationValues << new ConfigurationValue("controlSample", (bamControlMerged.fileStage as COFileStageSettings).sample.name)

            BaseFile sophiaControlFile = call(TOOL_SOPHIA, bamControlMerged, getISizesForBam(context, bamControlMerged, 0),"defaultReadLength=${controlDefaultReadLength}") as BaseFile
            BaseFile sophiaTumorFile = call(TOOL_SOPHIA, bamTumorMerged, getISizesForBam(context, bamTumorMerged, 1),"defaultReadLength=${tumorDefaultReadLength}") as BaseFile

            call("sophiaAnnotator", sophiaControlFile, sophiaTumorFile)
        } else { //NoControl!
            context.getConfiguration().getConfigurationValues().add(new ConfigurationValue(ANALYSIS_TAG, "tumorOnly"))
            BaseFile sophiaTumorFile = call(TOOL_SOPHIA, bamTumorMerged, getISizesForBam(context, bamTumorMerged, 0),"defaultReadLength=${tumorDefaultReadLength}") as BaseFile

            call("sophiaAnnotatorNoControl", sophiaTumorFile)
        }
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean initialBamCheck = super.checkExecutability(context)
        if (!initialBamCheck)
            return false

        BasicBamFile[] bamFiles = super.getInitialBamFiles(context)

        def configurationValues = context.getConfigurationValues()
        boolean isizesListIsSet = configurationValues.hasValue(ISIZES_LIST)
        boolean bamListIsSet = configurationValues.hasValue(BAMFILE_LIST)
        if(isizesListIsSet && ! bamListIsSet) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("Setting ${ISIZES_LIST} is only allowed if ${BAMFILE_LIST} is set as well."))
            return false
        }

        boolean isizesfound = true;
        for (int i = 0; i < bamFiles.size(); i++) {
            if (!bamFiles[i]) // Can be null for nocontrol
                continue

            def isizesFile = getISizesForBam(context, bamFiles[i], i).path
            boolean found = context.fileIsAccessible(isizesFile)
            if (!found) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("The insert sizes file ${isizesFile.path} could not be found."))
                isizesfound = false
            }
        }

        return isizesfound
    }
}
