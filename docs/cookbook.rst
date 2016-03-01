Cookbook
++++++++

This is a quick dive into how to build a Flake with this module.

* Create the Flake:

Import the growth module with ``import growth`` and create an instance of the
Flake: ``myFlake = growth.Flake()``

* Inspect the Flake:

  To view basic information, the representation of ``myFlake`` will display it's
  properties (**twinplanes**, **seed**, **iterations**, etc.). A visual
  representation is rendered in mayavi through ``myFlake.plot()``

* Grow atoms:

  Now we grow 20 atoms by calling the ``grow`` method on the Flake and passing
  the number as argument ``myFlake.grow(20)``. The growth method can be refined
  through various arguments, particularly the ``mode`` can be specified. As
  default we will grow with an exponential probability distribution according to
  our temperature in ``myFlake.temp``.


* Export results:

  After several iterations of inspecting and growing the Flake you can save the
  picture by passing ``save=True`` to the plot() method. This is then exported
  in output directory of the corresponding date.
  To save the coordinates in a `xyz`-file call the export method, optionally
  with a describing tag: ``myFlake.export('unexpected_phenomenon')``
  The resulting file is plaintext containing the coordinates and can further
  used or import into Blender for rendering.
