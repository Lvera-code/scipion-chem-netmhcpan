=========
CHANGES
=========

0.2.0
=====
- Merged scipion-chem-netmhciipan (MHC-II/T-helper, NetMHCIIpan-4.3) into
  this plugin as a second protocol, ``ProtNetMHCIIpanPromiscuity``: both
  tools share the same DTU academic-license installation shape, so there
  was no technical reason to keep them as two separate plugins. The two
  protocols themselves stay independent (see each docstring for why), and
  using only one does not require the other tool to be configured.

0.1.0
=====
- Initial release: NetMHCpan-4.2 MHC-I (HLA-A/B/C) cytotoxic T-cell
  promiscuity protocol, built as an independent plugin, mirroring the
  conventions established by scipion-chem-netmhciipan (MHC-II/T-helper),
  since NetMHCpan-4.2 also requires a separate academic-license
  installation (DTU Health Tech, not redistributable) but is a biologically
  distinct antigen presentation pathway that must never be merged with
  NetMHCIIpan's promiscuity figure.
