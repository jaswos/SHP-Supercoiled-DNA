from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules=[
    Extension("writhe_cython",
              sources=["writhe_cython.pyx"],
              libraries=["m"] # Unix-like specific
    )
]

setup(
  name = "writhe_cython",
  ext_modules = cythonize(ext_modules)
)
