from __future__ import print_function


class BrightliteWriter(object):

    def __init__(self, rc):
        self.rc = rc

    def write(self, lib, fname):
        header = ["TIME", "NEUT_PROD", "NEUT_DEST", "BUd"]
        lines = [h + " " + " ".join(map(str, lib[h])) for h in header]
        lines.extend([str(key).upper() + " " + " ".join(map(str, lib[key]))\
                      for key in lib.keys() if key not in header])

        with open(fname, "w") as f:
            f.write("\n".join(lines))
