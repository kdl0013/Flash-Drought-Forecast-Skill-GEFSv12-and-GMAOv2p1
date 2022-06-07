import distutils.core
import Cython.Build
import os, sys

cython_file = 'insert_cython_file'

#Language level 3 indicates python3
distutils.core.setup(
    ext_modules = Cython.Build.cythonize(f"{cython_file}"),include_dirs=os.path.abspath(os.path.dirname(sys.argv[0])),language_level=3)
