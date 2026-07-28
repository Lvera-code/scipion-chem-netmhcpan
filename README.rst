================================
NetMHCpan / NetMHCIIpan Scipion plugin
================================

Scipion framework plugin wrapping two DTU Health Tech tools for MHC
promiscuity prediction, merged 2026-07-28 (previously two separate
plugins, ``scipion-chem-netmhcpan`` and ``scipion-chem-netmhciipan``):
both wrap DTU academic-license tools with an identical installation
shape (manual download, edit the ``NMHOME`` line in the wrapper script,
point ``*_HOME`` at it in ``scipion.conf``), so there was no technical
reason to keep them as two plugins.

The plugin implements two protocols, kept separate on purpose (see each
protocol's own docstring): these are biologically distinct antigen
presentation pathways (professional antigen-presenting cell vs. any
nucleated cell; CD8+ cytotoxic vs. CD4+ helper; different binding-core
length rules and allele nomenclature), never merged into a single
promiscuity figure or a single form with an ``EnumParam`` toggle.

- ``ProtNetMHCpanPromiscuity`` (NetMHCpan-4.2): MHC-I (HLA-A/B/C)
  cytotoxic T-cell promiscuity, against a configurable panel of HLA-A/B/C
  alleles (23 alleles by default: 12 HLA-A/B supertypes of Sidney et al.
  2008 + 11 HLA-C alleles, chosen per Rasmussen et al. 2014 and IEDB's
  population-frequency criterion -- see ``constants.py`` for the full
  rationale and citations).
- ``ProtNetMHCIIpanPromiscuity`` (NetMHCIIpan-4.3): MHC-II (HLA-DR/DQ/DP)
  T-helper promiscuity, against a configurable panel of HLA-DR/DQ/DP
  alleles (the 27-allele IEDB reference panel by default).

Both take an arbitrary ``SetOfSequenceROIs`` (e.g. epitope candidates from
any upstream prediction/filtering protocol).

Neither tool is bundled with this plugin: both must be downloaded
separately (see below) and pointed to via ``scipion.conf``. Using only
ONE of the two protocols does not require the other tool to be set up.

================================
Download NetMHCpan-4.2 / NetMHCIIpan-4.3
================================

Both are **academic-use only** software (DTU Health Tech). Request them
from:

- NetMHCpan-4.2: https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/
- NetMHCIIpan-4.3: https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/

Each request form requires a valid **institutional/academic email
address** (a personal Gmail/Outlook address will be rejected) - use your
university or research center email. The download link is emailed to
that address.

For each tool, unzip it, then edit the ``NMHOME`` line inside its wrapper
script (``netMHCpan`` / ``netMHCIIpan``) with the absolute installation
path (a manual step required by DTU's own install instructions). Then, in
``scipion.conf``, set whichever of these you need:

.. code-block::

      NETMHCPAN_HOME = <path to the netMHCpan-4.2 folder>
      NETMHCIIPAN_HOME = <path to the netMHCIIpan-4.3 folder>

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-netmhcpan.git
            cd scipion-chem-netmhcpan
            scipion3 installp -p . --devel
