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

   To learn more about the run control file, see the `RunControl page <runcontrol.html>`_.
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
