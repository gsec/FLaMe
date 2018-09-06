This module contains defaults and global preference settings used across all simulations.

Before starting serious simulations you should set a main output path for the bulk growth
simulations, flake coordinate files and graphs.

It is easiest to set an environment variable and Flame outputs data relative to
that path.::

    $ export FLAME_OUTPUT=/your/path/here

Else the output path is set to ``../../../output`` relative to the ``settings.py`` file.

The options here can be seen as global settings, while in ``sim_params.yaml`` we define
simulations specific parameters or local settings, in each new simulation directory.

Distance => `DIFF_CAP`: This variable denotes the maximum distance (in atomic radii)
until which atoms are considered nearest neighbors when checked. For touching atoms this
is two atomic radii, while we give it here 2.1 to include small deviations due to
numerical rounding errors when calculating the distance (this is subject to further
considerations).
