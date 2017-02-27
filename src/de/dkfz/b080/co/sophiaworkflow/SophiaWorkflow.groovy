package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.b080.co.common.WorkflowUsingMergedBams
import de.dkfz.b080.co.files.BasicBamFile
import de.dkfz.roddy.config.ConfigurationValue
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

    private BaseFile getISizesForBam(BasicBamFile bam) {
        return BaseFile.constructManual(LibrariesFactory.getInstance().loadRealOrSyntheticClass("InsertSizesValueFile", BaseFile.class as Class<FileObject>), new BasicBamFile(bam)) as BaseFile
    }

    @Override
    boolean execute(ExecutionContext context, BasicBamFile bamControlMerged, BasicBamFile bamTumorMerged) {
        boolean fullWorkflow = bamControlMerged

        BaseFile insertSizesValueFileForTumor = getISizesForBam(bamTumorMerged)
        if (fullWorkflow) { // Control and tumor
            context.getConfiguration().getConfigurationValues().add(new ConfigurationValue("analysisTag", "tumorControl"))

            BaseFile insertSizesValueFileForControl = getISizesForBam(bamControlMerged)

            BaseFile sophiaControlFile = call("sophia", bamControlMerged, insertSizesValueFileForControl) as BaseFile
            BaseFile sophiaTumorFile = call("sophia", bamTumorMerged, insertSizesValueFileForTumor) as BaseFile

            call("sophiaAnnotator", sophiaControlFile, sophiaTumorFile)
        } else { //NoControl!
            context.getConfiguration().getConfigurationValues().add(new ConfigurationValue("analysisTag", "tumorOnly"))
            BaseFile sophiaTumorFile = call("sophia", bamTumorMerged, insertSizesValueFileForTumor) as BaseFile

            call("sophiaAnnotatorNoControl", sophiaTumorFile)
        }
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean initialBamCheck = super.checkExecutability(context)
        if(!initialBamCheck)
            return false

        BasicBamFile[] bamFiles = super.getInitialBamFiles(context)
        boolean isizesfound = true;
        for (int i = 0; i < bamFiles.size(); i++) {
            if(!bamFiles[i]) // Can be null for nocontrol
                continue

            def isizesFile = getISizesForBam(bamFiles[i]).path
            boolean found = context.fileIsAccessible(isizesFile)
            if(!found) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("The insert sizes file ${isizesFile.path} could not be found."))
                isizesfound = false
            }
        }

        return isizesfound
    }
}
