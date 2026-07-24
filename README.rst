================================
NetMHCpan Scipion plugin
================================

Scipion framework plugin wrapping NetMHCpan-4.2 (DTU Health Tech) for
MHC-I (HLA-A/B/C) cytotoxic T-cell promiscuity prediction. Deliberately
kept as an independent plugin from ``scipion-chem-netmhciipan`` (MHC-II/
T-helper): these are biologically distinct antigen presentation pathways
(professional antigen-presenting cell vs. any nucleated cell; CD8+
cytotoxic vs. CD4+ helper; different binding-core length rules), never
merged into a single promiscuity figure.

The plugin implements a single protocol, ``ProtNetMHCpanPromiscuity``,
which takes an arbitrary ``SetOfSequenceROIs`` (e.g. epitope candidates from
any upstream prediction/filtering protocol) and evaluates each peptide's
promiscuity across a configurable panel of HLA-A/B/C alleles (23 alleles by
default: 12 HLA-A/B supertypes of Sidney et al. 2008 + 11 HLA-C alleles,
chosen per Rasmussen et al. 2014 and IEDB's population-frequency criterion
-- see ``constants.py`` for the full rationale and citations).

NetMHCpan-4.2 is **not** bundled with this plugin: it must be downloaded
separately (see below) and pointed to via ``scipion.conf``.

================================
Download NetMHCpan-4.2
================================

NetMHCpan-4.2 is **academic-use only** software (DTU Health Tech). Request
it from:
https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/

The request form requires a valid **institutional/academic email address**
(a personal Gmail/Outlook address will be rejected) - use your university or
research center email. The download link is emailed to that address.

Unzip it, then edit the ``NMHOME`` line inside the ``netMHCpan`` wrapper
script with the absolute installation path (a manual step required by DTU's
own install instructions). Then, in ``scipion.conf``, set:

.. code-block::

      NETMHCPAN_HOME = <path to the netMHCpan-4.2 folder>

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-netmhcpan.git
            cd scipion-chem-netmhcpan
            scipion3 installp -p . --devel
