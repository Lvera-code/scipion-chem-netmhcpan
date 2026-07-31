# -*- coding: utf-8 -*-
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

# Merged 2026-07-28 (was two separate plugins, scipion-chem-netmhcpan and
# scipion-chem-netmhciipan): both wrap DTU Health Tech academic-license
# tools with an identical installation shape (manual download, edit the
# NMHOME line in the wrapper script, point *_HOME at it in scipion.conf) --
# there was no technical reason to keep them as separate plugins, only two
# separate protocol classes for two biologically distinct antigen
# presentation pathways (MHC-I cytotoxic vs MHC-II T-helper). See
# ProtNetMHCpanPromiscuity/ProtNetMHCIIpanPromiscuity docstrings for why
# they stay as two protocols, not one with an EnumParam.

NETMHCPAN_VERSION = '4.2'
NETMHCIIPAN_VERSION = '4.3'

# NetMHCpan-4.2 is academic-use only software (DTU Health Tech), not
# redistributable: it is never installed automatically. The user downloads
# it manually, edits the NMHOME line inside the netMHCpan wrapper script
# (a manual step required by DTU's own install instructions) and points
# NETMHCPAN_HOME to that folder in scipion.conf.
NETMHCPAN_DIC = {
    'name': 'NetMHCpan',
    'version': NETMHCPAN_VERSION,
    'home': 'NETMHCPAN_HOME',
    'binary': 'NETMHCPAN_BINARY_NAME',
}

DEFAULT_NETMHCPAN_BINARY_NAME = 'netMHCpan'

# NetMHCIIpan-4.3 is academic-use only software (DTU Health Tech), not
# redistributable: it is never installed automatically. The user downloads
# it manually, edits the NMHOME line inside the netMHCIIpan wrapper script
# (a manual step required by DTU's own install instructions) and points
# NETMHCIIPAN_HOME to that folder in scipion.conf.
NETMHCIIPAN_DIC = {
    'name': 'NetMHCIIpan',
    'version': NETMHCIIPAN_VERSION,
    'home': 'NETMHCIIPAN_HOME',
    'binary': 'NETMHCIIPAN_BINARY_NAME',
}

DEFAULT_NETMHCIIPAN_BINARY_NAME = 'netMHCIIpan'

# Both wrapper scripts ('netMHCpan'/'netMHCIIpan', the ones NETMHCPAN_HOME/
# NETMHCIIPAN_HOME point at) are '#!/bin/tcsh -f' scripts, not Python/bash --
# a real bug Blanca hit 2026-07-31: on a machine without tcsh installed
# system-wide, running them fails outright. tcsh itself (unlike the DTU
# binaries) is NOT license-restricted, so -- unlike NETMHCPAN_DIC/
# NETMHCIIPAN_DIC's binaries -- it IS safe to auto-install via conda.
# Shared by both tools (same generic dependency, no reason for two envs).
# The protocols invoke the wrapper scripts as 'tcsh <script> <args>'
# explicitly (see ProtNetMHCpanPromiscuity/ProtNetMHCIIpanPromiscuity's
# '_runMode') instead of executing them directly and relying on the
# shebang: the shebang hardcodes the absolute path '/bin/tcsh', which a
# conda-installed tcsh (living under the env's own prefix) would NOT
# satisfy even while the env is active -- only an explicit 'tcsh <script>'
# invocation actually resolves the interpreter via $PATH, which conda
# activation DOES affect.
TCSH_DIC = {
    'name': 'netmhc-tcsh',
    'version': '1.0',
    'home': 'NETMHC_TCSH_HOME',
    'activation': 'NETMHC_TCSH_ACTIVATION_CMD',
}

READ_URL = 'https://github.com/Lvera-code/scipion-chem-netmhcpan'
NETMHCPAN_DOWNLOAD_URL = 'https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/'
NETMHCIIPAN_DOWNLOAD_URL = 'https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/'

NETMHCPAN_NOINSTALL_WARNING = (
    'Installation could not be completed because the local NetMHCpan-4.2 '
    'installation has not been found. Due to academic license restrictions, '
    f'DTU Health Tech does not allow redistributing this package: download it '
    f'manually from {NETMHCPAN_DOWNLOAD_URL} (requires an academic account), '
    'edit the NMHOME line inside the netMHCpan wrapper script with the '
    'absolute installation path (a manual step required by DTU\'s own '
    'install instructions), and set NETMHCPAN_HOME in scipion.conf. Please '
    f'check the scipion-chem-netmhcpan README file for more details: {READ_URL}'
)

NETMHCIIPAN_NOINSTALL_WARNING = (
    'Installation could not be completed because the local NetMHCIIpan-4.3 '
    'installation has not been found. Due to academic license restrictions, '
    f'DTU Health Tech does not allow redistributing this package: download it '
    f'manually from {NETMHCIIPAN_DOWNLOAD_URL} (requires an academic account), '
    'edit the NMHOME line inside the netMHCIIpan wrapper script with the '
    'absolute installation path (a manual step required by DTU\'s own '
    'install instructions), and set NETMHCIIPAN_HOME in scipion.conf. Please '
    f'check the scipion-chem-netmhcpan README file for more details: {READ_URL}'
)

# Minimum footprint of an MHC-I binding peptide: below 8 aa there is no
# viable binding core for any known HLA-I allele. Kept as a hard floor, not
# a user parameter, since it reflects a binding-biology constraint rather
# than a tunable setting.
MIN_PEPTIDE_LENGTH_MHCI = 8

# Minimum footprint of the MHC-II binding core: NetMHCIIpan discards (or
# scores over a shorter, degraded core) peptides shorter than this. Kept as
# a hard floor, not a user parameter, since it reflects a binding-biology
# constraint rather than a tunable setting.
MIN_PEPTIDE_LENGTH_MHCII = 9

# Maximum safe length for the exact-peptide mode ('-p'): both the
# NetMHCpan-4.2 and NetMHCIIpan-4.3 binaries (Linux_x86_64) crash with a
# buffer overflow (SIGABRT) for longer inputs -- verified empirically
# against the reference panel: 55 aa OK, 57 aa crash for NetMHCpan (same
# buffer/architecture as NetMHCIIpan, whose limit falls in the same range;
# not a coincidence copied without checking, both share the same binary
# architecture). Longer sequences are routed to protein mode (sliding
# window) instead. A conservative safety margin (40, well below the
# discovered 55/57 aa boundary) is kept as a hard limit, not a user
# parameter, since it reflects a binary safety constraint rather than a
# tunable setting.
MAX_PEPTIDE_MODE_LENGTH = 40

# Candidate binding-core lengths used for NetMHCpan protein mode's internal
# sliding window ('-l' flag). MHC-I binding cores are canonically 8-11 aa
# (unlike NetMHCIIpan's fixed 15 aa window over a flexible 9 aa core), so
# this must be passed explicitly -- NetMHCpan does not infer it the way
# NetMHCIIpan's protein mode does.
PEPTIDE_LENGTHS = '8,9,10,11'

# Reference panel of 23 HLA-A/B/C alleles used to estimate broad population
# coverage when designing cytotoxic T-cell (MHC-I) epitopes. Never insert
# spaces between commas: NetMHCpan passes this string as-is to its '-a'
# allele parser, and a stray space breaks parsing of the following allele.
#
# HLA-A/B (12 alleles): representative of the 12 HLA class I supertypes of
# Sidney et al. 2008 ("HLA class I supertypes: a revised and updated
# classification", Immunome Res 4:2). That paper classifies ~700 HLA-A/-B
# alleles into 12 supertypes but does NOT cover HLA-C at all. The homologous
# MHC-I population-coverage panel most cited in the field (Weiskopf et al.
# 2013, PNAS, the "27-allele panel" popularized by IEDB for class I) is
# ALSO exclusively HLA-A/B (16 A + 11 B) -- so the historical absence of
# HLA-C from a panel like this one is a well-documented field-wide gap, not
# a deliberate omission.
#
# HLA-C (11 alleles): chosen for appearing consistently in two independent
# sources: (a) Rasmussen et al. 2014, J Immunol 193(10):4790-4801,
# "Uncovering the Peptide-Binding Specificities of HLA-C" (experimentally
# characterized, common-by-population-frequency HLA-C alleles); and (b)
# IEDB's >=1% global-population-frequency criterion, as applied by a 2020
# SARS-CoV-2 vaccine-design study that evaluated HLA-A/B/C jointly for
# maximum coverage (PMC7754929).
NETMHCPAN_REFERENCE_PANEL = (
    'HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,'
    'HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01,'
    'HLA-C01:02,HLA-C03:04,HLA-C04:01,HLA-C05:01,HLA-C06:02,HLA-C07:01,'
    'HLA-C07:02,HLA-C08:02,HLA-C12:03,HLA-C15:02,HLA-C16:01'
)

# IEDB reference panel of the 27 most representative HLA-DR/DQ/DP alleles,
# used to estimate broad population coverage when designing T-helper (MHC-II)
# epitopes. Never insert spaces between commas: NetMHCIIpan passes this
# string as-is to its '-a' allele parser, and a stray space breaks parsing
# of the following allele.
IEDB_REFERENCE_PANEL = (
    "DRB1_0101,DRB1_0301,DRB1_0401,DRB1_0405,DRB1_0701,DRB1_0802,DRB1_0901,"
    "DRB1_1101,DRB1_1201,DRB1_1302,DRB1_1501,DRB3_0101,DRB3_0202,DRB4_0101,DRB5_0101,"
    "HLA-DQA10501-DQB10201,HLA-DQA10501-DQB10301,HLA-DQA10301-DQB10302,"
    "HLA-DQA10401-DQB10402,HLA-DQA10101-DQB10501,HLA-DQA10102-DQB10602,"
    "HLA-DPA10201-DPB10101,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10401,"
    "HLA-DPA10301-DPB10402,HLA-DPA10201-DPB10501,HLA-DPA10201-DPB11401"
)
