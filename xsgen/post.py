"""The post-processing plugin for xsgen.
Post-processing Plugin API
===============
"""

import os
from xsgen.plugins import Plugin

class XSGenPlugin(Plugin):
    """This class cleans up after xsgen."""

    def execute(self, rc):
        """Leaves the reactor directory if --cwd was not specified."""

        if not rc.get('CWD'):
            os.chdir('..')
