#!/usr/bin/env python2
## encoding: utf-8
from mayavi import mlab
import growth

TL = 100
f = growth.Flake(-2, 2, trail=TL, temp=0)
coords  = f.plot(pipeline=True)

src = mlab.points3d(*coords, colormap="gist_ncar", scale_factor=0.1, vmin=0, vmax=15)
ms = src.mlab_source
fig = mlab.gcf()


@mlab.show
@mlab.animate(ui=False)
def anim():
  """Animate."""
  while 1:
    f.grow(TL)
    x, y, z, c = f.plot(pipeline=True)
    ms.reset(x=x, y=y, z=z, scalars=c)
    mlab.yaw(2)
    fig.scene.reset_zoom()
    yield

a = anim()
