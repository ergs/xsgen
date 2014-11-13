xsgen
=====
.. overview-start
Overview
--------

xsgen is a tool for computing multi-group neutron cross-section, burnup, and
multiplication factor (k\ :sub:`inf`\ ) as a function of reactor state, which
includes time/fluence, material properties, and reactor geometry.

.. image:: /xsgen-flow.svg

xsgen reads in reactor parameters, simulation parameters, and specific
quantities to keep track of from a run-control file. After validation, we use
these parameters to generate the set of possible reactor states. We group these
states into runs, by finding the ones that have the same initial conditions and
differ only by time.

For each run we have several timesteps, and for each timestep we run a neutron
transport code like OpenMC to find the multiplication factor k and the group
flux, φ\ :sub:`g`\ . We then feed these values to a transmutation code such as
ORIGEN2.2 to find the burnup and neutron production/destruction rates, as well
as the transmutation of the material itself. Once we have done this for all
timesteps, we write this out to libraries, theoretically of a variety of
formats.
.. overview-end

.. install-start

Installation
============

Dependencies
------------

xsgen depends on `PyNE <http://www.pyne.io>`_. You can find installation
instructions for PyNE `here <http://pyne.io/install.html>`_.

To install from source, simply download the source code from the
official GitHub repo and run ``setup.py``::

    git clone git://github.com/bright-dev/xsgen.git
    cd xsgen/
    python setup.py install --user

Currently, we require `OpenMC <http://mit-crpg.github.io/openmc/>`_
and `ORIGEN 2.2 <https://rsicc.ornl.gov/CustomerService.aspx>`_ for
burnup and criticality computations. Installation instructions, along
with a wealth of documentation, for OpenMC can be found `here
<http://mit-crpg.github.io/openmc/quickinstall.html>`_. ORIGEN 2.2 is
distributed by RSICC at the above link.

.. install-end
