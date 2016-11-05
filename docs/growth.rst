growth module
=============

The `growth` module is in charge to handle the attributes and the growing procedure of one
flake. It takes a grid object from the `grid` module and has methods to act on this grid.

The grid is accessed through the (i, j, k) cubic coordinates. Those can then be translated
through the grid.coord() method into actual cartesian (x, y, z) coordinates for the fcc 
structure including the twinplane errors.

The main idea for the growth module is that the lattice can be populated with `atoms` or 
`surface` sites. The wholeness of the atoms then concise the crystal we are growing, while
the surface is a layer of sites which may be populated in the next growth step.
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

# .. automodule:: growth
#     :members:
#     :show-inheritance:
