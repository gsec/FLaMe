****************
Modules Overview
****************

The :py:mod:`grid` module - a low level mapper
==============================================

.. include:: grid.rst

.. automodule:: flame.grid
    :members:

The :py:mod:`growth` module - manipulate the flake
==================================================

.. include:: growth.rst

.. automodule:: flame.growth
    :members:

The :py:mod:`settings` module - defaults and constants
======================================================

.. include:: settings.rst

.. automodule:: flame.settings
    :members:


Other modules
=============


flame.simulation module
-----------------------
Define the general structure of the simulations. We will create a hierarchy in the HDF
data format and group those elements into [simulation] > [TP state] > [samples]. A sample
is the smallest unit, a flake generated by the `growth` module.

.. automodule:: flame.simulation
    :members:

flame.paint module
------------------
Generate bokeh plots according to the given column and the flake data in HDF files.

.. automodule:: flame.paint
    :members:

flame.cli module
----------------
For quick usage instructions see :ref:`cookbook`.

.. automodule:: flame.cli
    :members:
