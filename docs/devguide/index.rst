Developer's Guide
=================

Overview of program flow
========================

First, ``xsgen`` reads in the initial conditions and other parameters
from the run control file, generating all possible states of the
permutation parameters.

Then, in ``xsgen.buk`` these states are sorted into runs, which share
all parameters except for burnup time. ``xsgen.buk`` then calls the
``generate_run`` method of the chosen solver for each run.

``generate_run`` returns a dictionary with the data we want to write
out to libraries. Depending on the solver its implementation may
differ, but this is how ``xsgen.openmc_origen`` works:

1. ``generate_run`` calls ``generate`` for each state in the run.
2. ``generate`` calls OpenMC for that reactor state, then reads the
   output to find the total neutron flux.
3. Using the neutron flux ``generate`` then calls
   ``run_all_the_origens``, running ORIGEN2.2 on the fuel
   material. This also runs ORIGEN once for each tracked nuclide - it
   starts with a material of 1000g of each nuclide. We then read the
   output to find the transmutation, burnup, and neutron
   production/destruction rates during this timestep.
4. ``generate`` then returns all this information for its timestep,
   and that gets merged into the previous timesteps by
   ``generate_run``.
5. Once all the timesteps are complete, ``generate_run`` returns the
   library data.

The library data is called ``libs`` in ``xsgen.buk``. It has significant nesting::
  libs
    |- "fuel"
    |- 1001
    |- 1003
    |- ... each 
  

**THIS IS WHAT LIBS LOOKS LIKE**


Plugins
=======

.. include:: ../plugins-module.rst

Writing Another Output Format
=============================

If everything you want is in libs, then just make a write() function
Otherwise you need to fix your solver to put more things in the libs.
