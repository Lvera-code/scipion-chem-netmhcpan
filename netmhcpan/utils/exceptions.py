"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class NetMHCpanExecutionError(Exception):
    """Failed to run NetMHCpan-4.2 locally: missing installation, failed/
    timed-out subprocess, or the output .xls was not generated / does not
    match the expected format."""
