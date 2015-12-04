import numpy as np
import pickle
# a = np.random.random((4, 4))
with open('grid.dat', 'rb') as f:
  a = pickle.load(f)
# a = np.asarray(a)
from mayavi.api import Engine
e = Engine()
e.start()
s = e.new_scene()
from mayavi.sources.api import ArraySource
src = ArraySource(scalar_data=a)
e.add_source(src)
from mayavi.filters.api import WarpScalar, PolyDataNormals
warp = WarpScalar()
e.add_filter(warp, obj=src)
normals = PolyDataNormals()
e.add_filter(normals, obj=warp)
from mayavi.modules.api import Surface
surf = Surface()
e.add_module(surf, obj=normals)
