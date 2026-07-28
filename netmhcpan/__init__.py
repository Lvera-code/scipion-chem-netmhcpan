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
This package contains protocols for MHC-I (HLA-A/B/C) cytotoxic T-cell
promiscuity prediction (local NetMHCpan-4.2) and MHC-II (HLA-DR/DQ/DP)
T-helper promiscuity prediction (local NetMHCIIpan-4.3). Merged into a
single plugin 2026-07-28 (previously scipion-chem-netmhcpan and
scipion-chem-netmhciipan): both wrap DTU Health Tech tools with an
identical installation shape, so there was no technical reason to keep
them as two plugins -- see ProtNetMHCpanPromiscuity/
ProtNetMHCIIpanPromiscuity for why they stay as two PROTOCOLS, not one
merged into a single EnumParam.
"""

import os

from pwchem import Plugin as pwchemPlugin

from .constants import (
    DEFAULT_NETMHCIIPAN_BINARY_NAME, DEFAULT_NETMHCPAN_BINARY_NAME, NETMHCIIPAN_DIC,
    NETMHCIIPAN_NOINSTALL_WARNING, NETMHCPAN_DIC, NETMHCPAN_NOINSTALL_WARNING,
)

_references = ['Reynisson2020', 'Nilsson2023']


class Plugin(pwchemPlugin):
    """NetMHCpan-4.2 and NetMHCIIpan-4.3 are both academic-use only
    software (DTU Health Tech): neither is ever installed automatically.
    See ``validateNetMHCpanInstallation``/``validateNetMHCIIpanInstallation``
    for what is checked per tool, and ``README.rst`` for the manual
    installation steps (identical shape for both: download, edit the
    NMHOME line inside the wrapper script, point *_HOME at it)."""

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(NETMHCPAN_DIC['home'], '')
        cls._defineVar(NETMHCPAN_DIC['binary'], DEFAULT_NETMHCPAN_BINARY_NAME)
        cls._defineVar(NETMHCIIPAN_DIC['home'], '')
        cls._defineVar(NETMHCIIPAN_DIC['binary'], DEFAULT_NETMHCIIPAN_BINARY_NAME)

    @classmethod
    def defineBinaries(cls, env):
        """No-op: neither NetMHCpan-4.2 nor NetMHCIIpan-4.3 is ever
        installed automatically (academic license, not redistributable).
        See ``validateNetMHCpanInstallation``/
        ``validateNetMHCIIpanInstallation``."""
        pass

    @classmethod
    def validateInstallation(cls):
        """Called by Scipion's plugin manager for the plugin as a whole:
        union of both tools' checks, informational only. Using only ONE of
        the two protocols does NOT require the other tool to be configured
        -- each protocol's own ``_validate()`` hooks into the
        tool-specific method instead (``validateNetMHCpanInstallation``/
        ``validateNetMHCIIpanInstallation``)."""
        return cls.validateNetMHCpanInstallation() + cls.validateNetMHCIIpanInstallation()

    @classmethod
    def validateNetMHCpanInstallation(cls):
        """Check that NetMHCpan-4.2 (MHC-I) is correctly installed. Returns
        a list of actionable error messages, empty if the installation is
        correct."""
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
            errors.append(NETMHCPAN_NOINSTALL_WARNING)
        return errors

    @classmethod
    def validateNetMHCIIpanInstallation(cls):
        """Check that NetMHCIIpan-4.3 (MHC-II) is correctly installed.
        Returns a list of actionable error messages, empty if the
        installation is correct."""
        errors = []

        home = cls.getNetMHCIIpanHome()
        if not home or not os.path.isdir(home):
            errors.append(f"NETMHCIIPAN_HOME is not set or does not exist: '{home}'.")
        else:
            binary = cls.getNetMHCIIpanBinaryPath()
            if not os.path.isfile(binary):
                errors.append(f"Could not find the local NetMHCIIpan-4.3 wrapper script at '{binary}'.")
            elif not os.access(binary, os.X_OK):
                errors.append(f"'{binary}' is not executable (chmod +x).")

        if errors:
            errors.append(NETMHCIIPAN_NOINSTALL_WARNING)
        return errors

    # ---------------------------------- Utils: NetMHCpan (MHC-I) -----------------

    @classmethod
    def getNetMHCpanHome(cls):
        return cls.getVar(NETMHCPAN_DIC['home'])

    @classmethod
    def getNetMHCpanBinaryPath(cls):
        home = cls.getNetMHCpanHome()
        if not home:
            return None
        return os.path.join(home, cls.getVar(NETMHCPAN_DIC['binary']))

    # ---------------------------------- Utils: NetMHCIIpan (MHC-II) --------------

    @classmethod
    def getNetMHCIIpanHome(cls):
        return cls.getVar(NETMHCIIPAN_DIC['home'])

    @classmethod
    def getNetMHCIIpanBinaryPath(cls):
        home = cls.getNetMHCIIpanHome()
        if not home:
            return None
        return os.path.join(home, cls.getVar(NETMHCIIPAN_DIC['binary']))
