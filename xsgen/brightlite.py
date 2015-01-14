from __future__ import print_function


class BrightliteWriter(object):

    def __init__(self, rc):
        self.rc = rc

    def write(self, libs, fname):
        rownames = ["TIME", "NEUT_PROD", "NEUT_DEST", "BUd"]
        for mat, matlib in libs.items():
            lines = [row + "   " + "   ".join(map(str, matlib[row]))
                     for row in rownames]
            nucs = matlib["tracked_nucs"]
            lines.extend([n + "   " + "   ".join(map(str, nucs[n]))
                          for n in nucs])
            with open(fname+"-"+str(mat), "w") as f:
                print("Writing out to {}...".format(fname+"-"+str(mat)))
                f.write("\n".join(lines))
