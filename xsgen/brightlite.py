from __future__ import print_function
import os
import shutil
from pyne import nucname
from math import pi
import numpy as np


class BrightliteWriter(object):

    def __init__(self, rc):
        self.rc = rc

    def write(self, libs, dirname):
        """Write out libraries to a directory.

        Parameters
        ----------
        libs : dict
            The reactor libraries gleaned from buk.
        dirname : str
            The output directory.
        """
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        rownames = ["TIME", "NEUT_PROD", "NEUT_DEST", "BUd"]
        for mat, matlib in libs.items():
            if isinstance(mat, int):
                fname = str(nucname.zzaaam(mat))
            else:
                fname = mat
            lines = [row + "   " + "   ".join(map(str, matlib[row]))
                     for row in rownames]
            nucs = matlib["tracked_nucs"]
            lines.extend(sorted([n + "   " + "   ".
                                 join(["{:.4g}".format(f) for f in nucs[n]])
                                 for n in nucs]))
            with open(os.path.join(dirname, fname + ".txt"), "w") as f:
                f.write("\n".join(lines))
        track_actinides = [n for n in nucs if nucname.znum(n) in nucname.act]
        with open(os.path.join(dirname, "manifest.txt"), "w") as f:
            f.write("\n".join([str(nucname.zzaaam(act)) for act in track_actinides]))
            f.write("\n")
        with open(os.path.join(dirname, "params.txt"), "w") as f:
            if self.rc.get("enrichment") is None:
                enrichment = self.rc.initial_heavy_metal.get(922350)
            else:
                enrichment = self.rc.enrichment
            if enrichment is not None:
                f.write("ENRICHMENT {}\n".format(enrichment))
            if self.rc.get("batches") is not None:
                f.write("BATCHES {}\n".format(self.rc.batches))
            if self.rc.get("pnl") is not None:
                f.write("PNL {}\n".format(self.rc.pnl))
            f.write("BURNUP {}\n".format(sum(libs["fuel"]["BUd"])))
            f.write("FLUX {:.0E}\n".format(np.mean(libs["fuel"]["phi_tot"][1:])))
        with open(os.path.join(dirname, "structural.txt"), "w") as f:
            clad_linear_density = pi * self.rc.clad_density * \
                (self.rc.clad_cell_radius ** 2 - self.rc.void_cell_radius ** 2)
            fuel_linear_density = pi * self.rc.fuel_density * \
                self.rc.fuel_cell_radius ** 2
            clad_frac = float(clad_linear_density / fuel_linear_density)
            cladrows = ["{} {:.8f}".format(nucname.zzaaam(n), f*clad_frac)
                        for n, f in self.rc.clad_material.comp.items()]
            f.write("\n".join(cladrows))
            f.write("\n")
        shutil.copyfile("TAPE9.INP", os.path.join(dirname, "TAPE9.INP"))
