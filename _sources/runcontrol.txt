Run Control files
=================

Much of the customization of an ``xsgen`` run is controlled through
the run control file. This is a Python file that consists of variable
assignments. For example, to set the ``reactor`` property of the run
control to ``"lwr1g"``, you would include this in the run control file::

  reactor = "lwr1g"

As a Python file, you can easily add your own comments and even set
parameters to Python expressions or objects from other imported
modules. A sample run control file can be found `here
<https://github.com/bright-dev/xslibs/blob/master/lwr1g.py>`_.

All plugins are able to update the command-line flags which set
options in the run control. You can see the individual `plugin pages
<plugins.html>`_ for details.

However, there are many options that are not set by command line -
these are reactor parameters that are loaded from the run control
file and validated by ``xsgen.pre``. To specify the run control file,
just add the flag ``--rc /path/to/file`` to your ``xsgen`` call.

If any parameters are missing, ``xsgen.pre`` will provide the following
defaults::

    defaultrc = {'formats': ('brightlite',),
                 'ui': False,
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
will certainly allow you to run xsgen.
