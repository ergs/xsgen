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

   ..
      To learn more about the run control file, see the `RunControl page<runcontrol.html>`_.
   For more information on command-line arguments, see the `plugins page <plugins.html>`_.


