from __future__ import print_function
import os
import shutil
from pyne import nucname
import numpy as np


class BrightliteWriter(object):

    def __init__(self, rc):
        self.rc = rc

    def write(self, libs, dirname):
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        rownames = ["TIME", "NEUT_PROD", "NEUT_DEST", "BUd"]
        for mat, matlib in libs.items():
            lines = [row + "   " + "   ".join(map(str, matlib[row]))
                     for row in rownames]
            nucs = matlib["tracked_nucs"]
            lines.extend(sorted([n + "   " + "   ".
                                 join(["{:.4g}".format(f) for f in nucs[n]])
                                 for n in nucs]))
            with open(os.path.join(dirname, str(mat)+".txt"), "w") as f:
                f.write("\n".join(lines))
        track_actinides = [n for n in nucs if nucname.znum(n) in nucname.act]
        with open(os.path.join(dirname, "manifest.txt"), "w") as f:
            f.write("\n".join([nucname.zzaaam(act) for act in track_actinides]))
        with open(os.path.join(dirname, "params.txt"), "w") as f:
            f.write("ENRICHMENT {}\n".format(self.rc.enrichment))
            f.write("BATCHES {}\n".format(self.rc.batches))
            f.write("PNL {}\n".format(self.rc.pnl))
            f.write("BURNUP {}\n".format(sum(libs["fuel"]["BUd"])))
            f.write("FLUX {:.0E}\n".format(np.mean(libs["fuel"]["phi_tot"][1:])))
        with open(os.path.join(dirname, "structural.txt"), "w") as f:
            coolrows = ("{} {:.8f}".format(nucname.zzaaam(n), f)
                        for n, f in self.rc.cool_material.comp.items())
            cladrows = ("{} {:.8f}".format(nucname.zzaaam(n), f)
                        for n, f in self.rc.clad_material.comp.items())
            f.write("\n".join(coolrows))
            f.write("\n")
            f.write("\n".join(cladrows))
        shutil.copyfile("TAPE9.INP", os.path.join(dirname, "TAPE9.INP"))

