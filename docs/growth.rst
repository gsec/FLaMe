This module is in charge to handle attributes and the growth procedure of a flake. It
instantiates a grid object from the :py:class:`Grid` class manages the atoms and surface
of the Flake. The ``Grid`` is used to translate the indexed atoms into real space with
proper coordinates and correct stacking permutation.

The grid is accessed through the (i, j, k) cubic coordinates. Those can then be translated
through the grid.coord() method into actual cartesian (x, y, z) coordinates for the fcc
structure including the twinplane errors.

To grow such a flake we have the :py:class:`Flake` class handling the atoms and surface
sites we use to grow the flake.
The lattice can be populated with `atoms` or `surface` sites. The wholeness of the atoms
then concise the crystal we are growing, while the surface is a layer of sites which may
be populated in the next growth step.

Here a simple (2D) ASCII illustration of the growth.

.. 

    * atom site    ``o`` 
    * surface site ``*`` 

::

        * ** *        * ** * *        * ** * *
        * oo *  ==>   * oo o *  ==>   * oo o *           etc..
        * ** *        * ** * *        * ** o *
                                         * * *

The growth can be done in three different modes:
  * random:
    Each surface site has the same probability to be populated in the next step.
  * deterministic:
    The sites with the most adjacent atoms will be populated next. If there are more then
    one with the most adjacent atoms, one is chosen randomly between those.
  * probabilistic:
    Each site is weighted through a specific function (e.g. exponential, Boltzmann, etc.)
    of the number of adjacent atoms and chosen according to that probability distribution.

The most interesting here, on which we also will focus in further discussions is the
probabilistic mode. The function we start with is an exponential governed by some kind of
'temperature'.
