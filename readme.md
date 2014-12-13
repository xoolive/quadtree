Quadtrees iterating on pairs of neighbouring items
==================================================

Implementation of a smart version of quadtrees specialised in tracking distance
between moving objects. When iterating over elements in the quadtree,
a decorator on the function to be iterated indicates whether the element might
have moved to a neighbouring subdivision.

# Installation

The easiest way to install the quadtree package is using pip
```
pip install git+git://github.com/xoolive/quadtree
```

# Usage

```python
import smartquadtree
```

You can refer to the `tutorial.iypnb` file for an introduction to the
possibilities of the package. If you are not familiar with the iPython
notebook, you can have a read-only access to the notebook through the
[viewer](http://nbviewer.ipython.org/github/xoolive/quadtree/blob/master/tutorial.ipynb).
