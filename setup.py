from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("smartquadtree",
              ["smartquadtree.pyx", "quadtree.cpp", "neighbour.cpp"],
              extra_compile_args=["-std=c++11"],
              language="c++")
]

setup(name="smartquadtree",
      version="0.1",
      author="Xavier Olive",
      author_email="xavier@xoolive.org",
      description="Implementation of quadtree for moving objects",
      ext_modules=cythonize(extensions),
      )
