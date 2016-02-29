Cookbook
++++++++

This is a quick dive into how to build a Flake with this module.

* Create the Flake:

Import the growth module with ``import growth`` and create an instance of the 
Flake: ``myFlake = growth.Flake()``

* Inspect the Flake:

To view basic information, calling the object ``myFlake`` will display it's
properties (**twinplanes**, **seed**, **iterations**, **etc.**).
A visual representation is rendered in mayavi through 
``myFlake.plot()``
