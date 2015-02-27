Run Control files
=================

Much of the customization of an ``xsgen`` run is controlled through
the run control file.  This is a Python file that consists of variable
assignments.  For example, to set the ``reactor`` property of the run
control to ``"lwr1g"``, you would include this in the run control file::

  reactor = "lwr1g"

As a Python file, you can easily add your own comments and even set
parameters to Python expressions or objects from other imported
modules.  A sample run control file can be found `here
<https://github.com/bright-dev/xslibs/blob/master/lwr1g.py>`_.

All plugins are able to update the command-line flags which set
options in the run control.  You can see the individual `plugin pages
<plugins.html>`_ for details.

However, there are many options that are not set by command line -
these are reactor parameters that are loaded from the run control
file and validated by ``xsgen.pre``.  To specify the run control file,
just add the flag ``--rc /path/to/file`` to your ``xsgen`` call.

If any parameters are missing, ``xsgen.pre`` will provide the following
defaults::

    defaultrc = {'formats': ('brightlite',),
                 'is_thermal': True,
                 'burn_times': [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
                 'group_structure': [10 ** x for x in range(1, -3, -1)],
                 'track_nucs': transmute,
                 'energy_grid': 'nuclide',
                 'fuel_density': 19.1,
                 'clad_density': 6.56,
                 'cool_density': 1.0,
                 'fuel_cell_radius': 0.7,
                 'void_cell_radius': 0.8,
                 'clad_cell_radius': 0.9,
                 'unit_cell_pitch': 1.5,
                 'unit_cell_height': 2,
                 'burn_regions': 1,
                 'fuel_specific_power': 1,
                 'fuel_material': {'U238': 0.96, 'U235': 0.04},
                 'solver': "openmc+origen",
                 'reactor': "lwr",
                 'k_cycles': 20,
                 'k_cycles_skip': 10,
                 'k_particles': 1000,
                 }

There is no guarantee that these are particularly physical, but they
will certainly allow you to run xsgen. Here is what they all mean:

* ``formats`` is a tuple containing the desired output
  formats. Currently, only ``brightlite`` is available.
* ``is_thermal`` is used to determine whether we can use EAF data, and
  which ORIGEN call to make. When ``True``, our reactor is
  thermal. When ``False``, it's fast.
* ``burn_times`` is an array of the burn times that we want to get
  data for. The units are days.
* ``group_structure`` is an array of the energy bins you want to use,
  in MeV.
* ``track_nucs`` is an array of the nuclides you would like to
  generate libraries for. Usually it will be something like
  [922350000, 922380000].
* ``energy_grid`` can be either "nuclide" or "unionized", which gets
  written to OpenMC's settings. Beware - unionized energy grid can
  cause OpenMC to run out of memory. If OpenMC is crashing with exit
  code -9, you might want to make sure ``energy_grid`` is "nuclide".
* ``fuel_density``, ``clad_density``, and ``cool_density`` have units
  of [g/cm^3].
* ``fuel_cell_radius``, ``void_cell_radius``, ``clad_cell_radius``,
  ``unit_cell_pitch``, and ``unit_cell_height`` have units of [cm].
* ``burn_regions`` is the number of annular burn regions.
* ``fuel_specific_power`` has units of [W/g].
* ``fuel_material`` is a dictionary with the mass composition of the
  fuel material.
* ``solver`` is the transport and transmutation engine - we only
  support "openmc+origen" right now.
* ``reactor`` is the name of the reactor. It's used for the default
  output directory.
* ``k_cycles`` is the number of cycles to run the transport for.
* ``k_cycles_skip`` is the number of cycles to skip initially.
* ``k_particles`` is the number of particles to run in each cycle.

Additional parameters with no defaults include the following. To
specify them you can put them in the run control file.

* ``plugins`` is a list of plugins you want to run.
* ``burn_time`` is the maximum burn time you want to calculate
  libraries for. You can set this, together with ``time_step``,
  instead of specifying each burn time individually in ``burn_times``.
  However ``burn_times`` will overrride this. This has units of [days].
* ``time_step`` is the time step, in units of [days], if using
  ``burn_time``. The burn times will end up counting up from 0 to
  ``burn_time`` in steps of ``time_step`` days.
* ``initial_heavy_metal`` is the mass fraction distribution of the initial heavy metal.
* ``fuel_chemical_form`` is the atomic fraction distribution of the
  initial fuel loading. One of the keys should be "IHM" for the
  initial heavy metal.
* ``temperature`` has units of [K] and is the cross-section
  temperature we use. It should be a multiple of 300K.

The next few are for the Bright-lite input files.

* ``enrichment`` is the enrichment for the ``params.txt`` Bright-lite
  input file.  If it's not specified, we try to find the proportion of
  922350 in the initial heavy metal.  If that's not there, then we
  don't write out enrichment at all.
* ``batches`` is the number of batches to put in the ``params.txt``
  Bright-lite input file.  If not specified, then we won't write
  anything about batches to ``params.txt``.
* ``pnl`` The PNL to put in the ``params.txt`` Bright-lite input file.
  If not specified, then we won't write anything about PNL to
  ``params.txt``.

And this last one is for configuring OpenMC.

* ``openmc_cross_sections`` is the path to your OpenMC cross-sections.
  If not specified, we attempt to read the $CROSS_SECTIONS environment
  variable.
