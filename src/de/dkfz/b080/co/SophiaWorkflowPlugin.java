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

package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 * * TODO Recreate class. Put in dependencies to other workflows, descriptions, capabilities (like ui settings, components) etc.
 */
public class SophiaWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "2.2.0";
    public static final String CURRENT_VERSION_BUILD_DATE = "Tue Feb 26 12:50:40 CET 2019";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}

