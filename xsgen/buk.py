"""Plugin that runs burnup-criticality calculation.
"""
from __future__ import print_function
import os

from xsgen.plugins import Plugin
from xsgen.utils import RunControl, NotSpecified
from xsgen.openmc_origen import OpenMCOrigen

SOLVER_ENGINES = {'openmc+origen': OpenMCOrigen}

class XSGenPlugin(Plugin):

    requires = ('xsgen.pre',)

    defaultrc = RunControl(
        solver=NotSpecified,
        openmc_cross_sections=NotSpecified,
        )

    rcdocs = {
        'openmc_cross_sections': 'Path to the cross_sections.xml file for OpenMC',
        'solver': ('The physics codes that are used to solve the '
                   'burnup-criticality problem and compute cross sections and '
                   'transmutation matrices.'),
        }

    def update_argparser(self, parser):
        parser.add_argument('--solver', dest='solver', help=self.rcdocs['solver'])
        parser.add_argument("--openmc-cross-sections", dest="openmc_cross_sections",
            help=self.rcdocs['openmc_cross_sections'])

    def setup(self, rc):
        self._ensure_omcxs(rc)

        # do after all other values have been setup
        if rc.solver is NotSpecified:
            raise ValueError('a solver type must be specified')
        rc.engine = SOLVER_ENGINES[rc.solver](rc)

    def same_except_burnup_time(self, state1, state2):
        if len(state1) != len(state2):
            raise ValueError("States have unequal number of perturbation paramaters.")
        for index in range(len(state1)):
            if state1._fields[index] == 'burn_times':
                continue
            if state1[index] != state2[index]:
                return False
        return True

    def execute(self, rc):
        runs = []
        for state in rc.states:
            already_existed = False
            for run in runs:
                if self.same_except_burnup_time(run[0], state):
                    run.append(state)
                    already_existed = True
            if not already_existed:
                runs.append([state])
        rc.runs = runs

        for run in rc.runs:
            lib = rc.engine.generate_run(run)
            for writer in rc.writers:
                writer.write(lib)

    #
    # ensure functions
    #

    def _ensure_omcxs(self, rc):
        if rc.openmc_cross_sections is not NotSpecified: # which means Specified
            rc.openmc_cross_sections = os.path.abspath(rc.openmc_cross_sections)
        elif 'CROSS_SECTIONS' in os.environ:
            rc.openmc_cross_sections = os.path.abspath(os.environ['CROSS_SECTIONS'])
        else:
            rc.openmc_cross_sections = None
            
