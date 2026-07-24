# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains a protocol for MHC-I (HLA-A/B/C) cytotoxic T-cell
promiscuity prediction using a local NetMHCpan-4.2 installation.
"""

import os

from pwchem import Plugin as pwchemPlugin

from .constants import DEFAULT_BINARY_NAME, NETMHCPAN_DIC, NOINSTALL_WARNING

_references = ['Reynisson2020']


class Plugin(pwchemPlugin):
    """NetMHCpan-4.2 is academic-use only software (DTU Health Tech): it is
    never installed automatically. See ``validateInstallation`` for what is
    checked and ``README.rst`` for the manual installation steps."""

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(NETMHCPAN_DIC['home'], '')
        cls._defineVar(NETMHCPAN_DIC['binary'], DEFAULT_BINARY_NAME)

    @classmethod
    def defineBinaries(cls, env):
        """No-op: NetMHCpan-4.2 is never installed automatically (academic
        license, not redistributable). See ``validateInstallation``."""
        pass

    @classmethod
    def validateInstallation(cls):
        """Check that this plugin's requirements are met. Called by Scipion's
        plugin manager. Returns a list of actionable error messages, empty if
        the installation is correct."""
        errors = []

        home = cls.getNetMHCpanHome()
        if not home or not os.path.isdir(home):
            errors.append(f"NETMHCPAN_HOME is not set or does not exist: '{home}'.")
        else:
            binary = cls.getNetMHCpanBinaryPath()
            if not os.path.isfile(binary):
                errors.append(f"Could not find the local NetMHCpan-4.2 wrapper script at '{binary}'.")
            elif not os.access(binary, os.X_OK):
                errors.append(f"'{binary}' is not executable (chmod +x).")

        if errors:
            errors.append(NOINSTALL_WARNING)
        return errors

    # ---------------------------------- Utils -----------------------------------

    @classmethod
    def getNetMHCpanHome(cls):
        return cls.getVar(NETMHCPAN_DIC['home'])

    @classmethod
    def getNetMHCpanBinaryPath(cls):
        home = cls.getNetMHCpanHome()
        if not home:
            return None
        return os.path.join(home, cls.getVar(NETMHCPAN_DIC['binary']))
