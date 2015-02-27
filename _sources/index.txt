.. xsgen documentation master file, created by
   sphinx-quickstart on Fri Oct 31 11:35:50 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

XSGen: Cross Section Generator for Bright-lite Reactors
=======================================================
xsgen is a tool for computing multi-group neutron cross-section, burnup, and
multiplication factor (k\ :sub:`inf`\ ) as a function of reactor state, which
includes time/fluence, material properties, and reactor geometry.

We use a plugin-based system to allow a flexible approach with respect to
underlying physics engines.

The source code can be found at the `GitHub project site`_.

For a quick install from source, please clone from the official repo::

    git clone git://github.com/bright-dev/xsgen.git
    cd xsgen/
    python setup.py install --user

Now you should be able to run it simply with::

    xsgen

Release Notes
---------------
**0.1: 2/27/2015**
  * We can now run OpenMC and Origen together to generate Bright-lite
    input files!
  Developers:

  * Anthony Scopatz
  * John Xia

Contents:
---------

**Usage:**

.. toctree::
   :maxdepth: 1

   overview
   install
   usersguide
   runcontrol
   plugins
   api/index
   devguide/index

..
   Indices and tables
   ------------------

   * :ref:`genindex`
   * :ref:`modindex`
* :ref:`search`

.. _GitHub project site: https://github.com/bright-dev/xsgen
