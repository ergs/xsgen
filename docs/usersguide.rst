xsgen User's Guide
==================

1. Installation
===============

   .. include:: ../README.rst
      :start-after: .. install-start
      :end-before: .. install-end

   For this example, we will be using the OpenMC/ORIGEN physics code.


2. Setting Up the Run Control File
==================================

   For this example, we'll use the `xslibs
   <https://github.com/bright-dev/xslibs>`_ project repository by
   Anthony Scopatz. This provides some useful run-control files for
   light water reactors. ::

     git clone git@github.com:bright-dev/xslibs.git
     cd xslibs
     xsgen --rc lwr1g.py

   Let's look at what's inside ::

     reactor = "lwr1g"

   This is what you want to call the reactor. Here we call it
   ``lwr1g``. The output files will all go into a ``build-lwr1g``
   directory unless you specify otherwise with the ``--outdirs`` flag.

   ``xsgen`` uses a plugin-based system to run various chunks of
   code. ``xsgen.base`` is always included.  Next we see a list of the
   plugins we are using::

     plugins = ['xsgen.pre', 'xsgen.buk']

   ``xsgen.pre`` validates the run control file and adds in sane
   defaults for anything we are missing.  ``xsgen.buk`` takes care of
   the burnup and criticality calculations.  Additional plugins can be
   specified with the ``--plugins`` command-line flag.  There are a
   few more parameters specifying different aspects of the burnup and
   criticality calculations::

     solver = 'openmc+origen'
     formats = ('brightlite,')
     burn_regions = 1     # Number of burnup annular regions.
     burn_time = 365*10  # Number of days to burn the material [days]
     time_step = 100      # Time step by which to increment the burn [days]
     burn_times = [0, 3]
     burn_times.extend(range(100, 4001, 100))  # we now have [0, 3, 100, 200 .. 4000]

   Here, ``solver`` refers to the engine we use for the transmutation
   and transport. We currently only have one where we use OpenMC for
   neutron transport and Origen 2.2 for transmutation, named
   "openmc+origen." ``formats`` is an iterable with the various output
   formats you intend to generate libaries for. We currently only
   support `Bright-lite
   <https://github.com/FlanFlanagan/Bright-lite/tree/timestep/>format`_.
   ``burn_regions`` is, as the comment says, the number of burnup
   annular regions.

   You can specify your burn times in two different ways. You can
   either specify a ``burn_time`` and a ``time_step``, which will make
   the burn times start at zero and go to ``burn_time`` in increments
   of ``time_step``. Or you can set ``burn_times`` to a list of times.
   The burn times are in units of [days]. ``burn_times`` will take
   precedence over ``burn_time`` and ``time_step``. Additionally, if
   no burn times are specified, we default to ``[0, 100, ... , 900]``.

   Note: if you intend to use these libraries for Bright-lite, you
   should have a burn time at 3 days. This gets us a non-zero but
   small fluence, which helps with interpolation. You should also go
   out to at least 7000 days, otherwise you will not get enough
   fluence.

   Next we have settings for the geometry and fuel::

     fuel_cell_radius = 0.410
     void_cell_radius = 0.4185
     clad_cell_radius = 0.475
     unit_cell_pitch  = 0.65635 * 2.0
     unit_cell_height = 10.0

     fuel_density = 10.7
     clad_density = 5.87                         # Cladding Density
     cool_density = 0.73                         # Coolant Density

     fuel_specific_power = 40.0 / 1000.0   # Power garnered from fuel [W / g]

   The distance measurements here are all in [cm].  The densities are in
   [g/cm^3].  The fuel-specific power is in [W/g]. ::

     # LEU
     initial_heavy_metal = {     # Initial heavy metal mass fraction distribution
     922350: 0.04,
     922380: 0.96,
     }

     # UOX
     fuel_chemical_form = {                 # Dictionary of initial fuel loading.
     80160: 2.0,
     "IHM": 1.0,
     }

   Here we tell ``xsgen`` that our initial heavy metal is 4% U235 and
   96% U238.  ``fuel_chemical_form`` sets the initial fuel - here we
   see that we have one atom of IHM to every two atoms of oxygen, or
   UO_2.

   Next we set a few OpenMC parameters::

     k_particles   = 500       # Number of particles to run per kcode cycle
     k_cycles      = 130       # Number of kcode cycles to run
     k_cycles_skip = 30        # Number of kcode cycles to run but not tally at the begining.

   These set the various parameters relating to the transport - how
   many particles in each cycle, how many cyclues to run, and how many
   initial cycles to skip.

   For more information on the command-line arguments that various
   plugins supply, see the `plugins page <plugins.html>`_.

   Now ``xsgen`` will look at your run control, figure out how many
   different sets of initial conditions it needs to take into account,
   and calculate a separate run with OpenMC and Origen for each set of
   initial conditions. Each run will be output to a different
   directory.

   ``xsgen`` will build everything by default in the ``build-lwr1g``
   directory. In this directory will be separate directories which are
   named by hashes of each state. Inside these, there are ``omc``
   directories and ``origen<nuclide id>`` directories, which hold the
   input and output files for OpenMC and ORIGEN runs. The Bright-lite
   output files will be in ``brightliteN`` directories, where N is the
   run id.

   This is what the directory structure should look like::

     xslibs/
       |- build-lwr1g
            |- -4743197135481227544
                 |- omc
                      |- OpenMC input files go here
                 |- origen10010000
                      |- ORIGEN input files go here
                 |- origen10010000
                      |- ORIGEN input files go here
                 |- ...
                 |- origenfuel
                      |- ORIGEN input files for full fuel go here
            |- -4743197135546938812
            |- ...
            |- brightlite0
                 |- Bright-lite output files go here


3. Using the Libraries With Bright-lite
=======================================

   First we must install Bright-lite. Clone the repository
   `here<https://github.com/FlanFlanagan/Bright-lite/tree/timestep/>`_,
   then run ``python install.py``.

   To use the libraries with Bright-lite, we need to copy the folder
   ``brightlite0`` to
   ``${CYCLUS_INSTALL_DIR}/share/brightlite/brightlite0``.  Then we
   need to copy ``eg01-eg23.json`` from the ``xslibs`` repository into
   that folder. We also need to copy ``input/FR_reprocess.xml`` into
   ``.../brightlite/brightlite0/hist``. Now we have all the files in
   the right place and can run Cyclus::

     $ cyclus eg01-eg23.json

     New brightlite0 reactor starting up - target burnup = 45
     /home/john/cyclus/install/share/brightlite/brightlite0/manifest.txt
     [...]
     CR: 191.444
     Cycle Months:  10
     Cycle Days: 10.9494
     BURNUP: 45.1533

     Status: Cyclus run successful!
     Output location: cyclus.sqlite
     Simulation ID: de741da0-39d3-4933-a6c5-0f15e78813e5
