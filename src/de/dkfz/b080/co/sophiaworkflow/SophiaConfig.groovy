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

import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.core.ExecutionContext
import groovy.transform.CompileStatic

@CompileStatic
class SophiaConfig {

    private ExecutionContext context
    private Configuration configuration

    SophiaConfig(ExecutionContext context) {
        this.context = context
        this.configuration = context.getConfiguration()
    }
}
