from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("smartquadtree",
              ["smartquadtree.pyx", "quadtree.cpp", "neighbour.cpp"],
              extra_compile_args=["-std=c++11"],
              language="c++")
]

def get_long_description():
    import codecs
    with codecs.open('tutorial.rst', encoding='utf-8') as f:
        readme = f.read()
    return readme

setup(name="smartquadtree",
      version="1.0",
      author="Xavier Olive",
      author_email="xavier@xoolive.org",
      description="Implementation of quadtrees for moving objects",
      long_description=get_long_description(),
      license="MIT",
      url="https://github.com/xoolive/quadtree",
      ext_modules=cythonize(extensions),
      )

# Producing long description
# ipython nbconvert tutorial.ipynb --to rst
# Then manually edit paths to images to point to github
