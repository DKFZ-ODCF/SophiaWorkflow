package de.dkfz.b080.co.sophiaworkflow

import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.core.ExecutionContext
import groovy.transform.CompileStatic

/**
 * Created by heinold on 14.02.17.
 */
@CompileStatic
class SophiaConfig {

    private ExecutionContext context
    private Configuration configuration

    SophiaConfig(ExecutionContext context) {
        this.context = context
        this.configuration = context.getConfiguration()
    }
}
