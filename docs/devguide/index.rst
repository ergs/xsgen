Developer's Guide
=================

Overview of program flow
------------------------

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
          |- "material" : a list of pyne.material.Material's per timestep
          |- "tracked_nucs"
              |- <nuclide>: a list of <nuclide>'s mass fractions per timestep
              |- <nuclide>: a list of <nuclide>'s mass fractions per timestep
              |- ...
          |- "TIME": a list of timesteps, e.g [0, 3, 100, 200, ...]
          |- "NEUT_DEST": a list of neutron destruction rates per timestep
          |- "NEUT_PROD": a list of neutron production rates per timestep
          |- "BUd": a list of burnup per timestep
          |- "phi_tot": a list of neutron flux per timestep
      |- 1001
      |- 1003
      |- ... each nuclide we track has its own dictionary with the above structure


Writing To Another Output Format
--------------------------------

One of the things you may want to do is output to other formats.

To do this, you need to make ``some_format.py``, which will contain a
class with a ``write`` method. The ``write`` method will take ``libs``
and ``dirname`` as additional parameters.

The ``libs`` described above will be passed in, so you just need to
turn that dictionary into whatever output format you would like. If
your desired format needs information not in ``libs`` you will need to
do some additional work to make the solver put the right data in
``libs``.

Plugins
-------

.. include:: ../plugins-module.rst
