/*
 * Copyright (c) 2017 eilslabs.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/COWorkflowsBasePlugin/LICENSE.txt).
 */
package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.roddy.execution.io.BaseMetadataTable
import groovy.transform.CompileStatic

@CompileStatic
/** Metadata table to guide Sophia analyses.
 *
 * fileCol: BAM file path
 * datasetCol: Name of the dataset, i.e. PID
 * sampleCol: Name of the sample
 * sampleClassCol: Type fo the sample; tumor/probe vs. control/reference
 * defaultReadLengthCol: Read length.
 * medianIsizeCol: Median insert size for the BAM file.
 * stdIsizePercentageCol: ???
 * properPairRatioCol: Ratio of properly paired alignment pairs / total paired alignments.
 *
 */
class SophiaMetadataTable extends BaseMetadataTable {

    enum SampleClass {
        Reference,
        Probe

        static SampleClass fromString(String name) {
            switch (name.toLowerCase()) {
                case "reference":
                    Reference
                    break
                case "control":
                    Reference
                    break
                case "probe":
                    Probe
                    break
                case "tumor":
                    Probe
                    break
                default:
                    throw new IllegalArgumentException("Unallowed sample class '${name}'. Use reference/control or probe/tumor.")
            }
        }

        @Override
        String toString() {
            if (this == Reference)
                return "control"
            else (this == Probe)
                return "tumor"
        }
    }

    SophiaMetadataTable(BaseMetadataTable baseMetadataTable) {
        super(baseMetadataTable)
    }

    static void ensureDefault(Map<String, String> idmap, String key, String defaultValue) {
        if (!idmap.containsKey(key))
            idmap[key] = defaultValue
    }

    private void assertUniqueFile() {
        Map<String, Integer> tooFrequentFiles = records.countBy {
            it.get(INPUT_TABLE_FILE)
        }.findAll { file, count ->
            count > 1
        }
        if (tooFrequentFiles.size() > 0) {
            throw new RuntimeException("Files occur too often in input table: ${tooFrequentFiles}")
        }
    }

    @Override
    BaseMetadataTable subsetByColumn(String columnName, String value) {
        return new SophiaMetadataTable(super.subsetByColumn(columnName, value))
    }

}
