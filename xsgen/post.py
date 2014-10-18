import os
from xsgen.plugins import Plugin

class XSGenPlugin(Plugin):
    def execute(self, rc):
        # Clean up
        if not rc.get('CWD'):
            os.chdir('..')
