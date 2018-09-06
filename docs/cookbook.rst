.. _cookbook:

********
Cookbook
********


Installation
------------
  Install under Python 3. Please note that under Python 3 it is no longer possible to
  create a 3D view of the flake due to the lack of the ``mayavi`` package.
  Switch to main directory and..

  ... first install requirements...

  ``$ pip install -r requirements.txt``

  ... and then the flame module:

  ``$ python setup.py install``


Single Flake
------------
This is run in a python console or as script.

* Create the Flake:
  Import the growth module e.g. ``from flame import growth`` and create an instance of the
  Flake: ``myFlake = growth.Flake()``

* Inspect the Flake:
  To view basic information, the representation of ``myFlake`` will display it's
  properties (twinplanes, seed, iterations, etc.). A visual representation is
  rendered in mayavi through ``myFlake.plot()``

* Grow atoms:
  Now we grow 20 atoms by calling the ``grow`` method on the Flake and passing the number
  as argument ``myFlake.grow(20)``. The growth method can be refined through various
  arguments, particularly the ``mode`` can be specified. As default we will grow with an
  exponential probability distribution according to our temperature in ``myFlake.temp``.

* Export results(**OUTDATED**):
  After several iterations of inspecting and growing the Flake you can save the picture by
  passing ``save=True`` to the plot() method. This is then exported in output directory of
  the corresponding date. To save the coordinates in a `xyz`-file call the export method,
  with a describing tag: ``myFlake.export_coordinates('unexpected_phenomenon')`` The
  resulting file is plaintext containing the coordinates and can be used to import the
  coordinates into Blender for rendering or further use.

Automated Flake generation
--------------------------
This is run in a shell calling the program with arguments.

* Create a simulation:
  - create new directory with the desired name
  - change into the directory
  - run ``$ flame create``
  This will create a new file name ``sim_params.yaml`` containing conservative default
  settings for a new simulation. Those can be edited as needed. The most relevant is the
  function and the values for creating twin planes. This is used to sweep this parameter
  as desired.

* Run simulation:
  After changing the settings to your desire run ``$ flame run`` and wait for the
  simulation to run. A HDF5 file is created with the name of the project and a random salt
  to make it possible to run multiple simulations with the same parameters. The simulation
  grows flakes with the given parameters and at each growth step defined by ``snapshot
  interval`` the output of the ``geometry()`` method is then saved.

* View results:
  After the simulation is done, run ``$ flame paint <column1> <column2> <...>`` in the
  project folder. Where the columns are the aspects of interest of the ``geometry()``
  method. Each result is then rendered and saved to a file named ``<column>.html`` and can
  be viewed, scrolled and panned in a regular browser.
