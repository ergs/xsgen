xsgen User's Guide
==================

1. Installation

   .. include:: ../README.rst
      :start-after: .. install-start
      :end-before: .. install-end

   For this example, we will be using the OpenMC/ORIGEN physics code.


2. Usage

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
   code. Next we see a list of the plugins we are using::

     plugins = ['xsgen.pre', 'xsgen.buk']

   ``xsgen.pre`` validates the run control file and adds in sane
   defaults for anything we are missing. ``xsgen.buk`` takes care of
   the burnup and criticality calculations. Additional plugins can be
   specified with the ``--plugins`` command-line flag . There are a
   few more parameters specifying different aspects of the burnup and
   criticality calculations::

     solver = 'openmc+origen'
     formats = ('brightlite,')
     burn_regions = 1     # Number of burnup annular regions.
     burn_time = 365*10  # Number of days to burn the material [days]
     time_step = 100      # Time step by which to increment the burn [days]
     batches = 3

   Here, ``solver`` refers to the engine we use for the transmutation
   and transport. We currently only have one where we use OpenMC for
   neutron transport and Origen 2.2 for transmutation, named "openmc+origen."

   For more information on command-line arguments, see the `plugins page <plugins.html>`_.

   ``xsgen`` will build everything by default in the ``build-lwr1g``
   directory. In this directory will be separate directories which are
   named by hashes of each state. Inside these, there are ``omc``
   directories and ``origen<nuclide id>`` directories, which hold the
   input and output files for OpenMC and ORIGEN runs.

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
