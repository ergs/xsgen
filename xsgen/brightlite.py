from __future__ import print_function

class BrightliteWriter(object):

    def __init__(self, rc):
        self.rc = rc

    def write(self, state, lib):
        print('lib is {}'.format(lib))
