/*
 * Copyright (C) 2018 Michael Heinold, Philip R. Kensche and DKFZ Heidelberg
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.b080.co.common.WorkflowUsingMergedBams
import de.dkfz.b080.co.files.BasicBamFile
import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationError
import de.dkfz.roddy.config.ConfigurationValue
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues as RecCont
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.knowledge.files.BaseFile
import groovy.transform.CompileStatic

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
            if (cvalue.empty || cvalue == "")
                throw new ConfigurationError("Configuration value '${cvalueName}' is unset", config)
            new MapEntry(varname, cvalue)
        }
        result
    }

    @Override
    boolean execute(ExecutionContext context, BasicBamFile bamControlMerged, BasicBamFile bamTumorMerged) {
        boolean fullWorkflow = bamControlMerged

        Configuration configuration = context.configuration
        RecCont configurationValues = configuration.configurationValues

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

    @Override
    boolean checkExecutability() {
        boolean initialBamCheck = super.checkExecutability()
        if (!initialBamCheck)
            return false

        loadInitialBamFilesForDataset()

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
