# distutils: language = c++
# -*- coding: utf-8 -*-

""" An implementation of quadtrees -- iterating on pairs of neighbouring items

A quadtree is a tree data structure in which each node has exactly four
children. It is a particularly efficient way to store elements when you need
to quickly find them according to their x-y coordinates.

A common problem with elements in quadtrees is to detect pairs of elements
which are closer than a definite threshold.

The proposed implementation efficiently addresses this problem.

"""

from libcpp cimport bool as cppbool
from libcpp.vector cimport vector
from cpython.ref cimport Py_INCREF, Py_DECREF
from cython.operator cimport dereference as deref, preincrement as inc


cdef extern from *:
    cdef cppclass cpp_vector "std::vector" [T]:
        cppclass const_iterator:
            T& operator*()
            const_iterator operator++()
            bint operator!=(const_iterator)
        const_iterator begin()
        const_iterator end()

cdef extern from "quadtree.h":
    cdef cppclass PolygonMask:
        PolygonMask(vector[double], vector[double], int)
    cdef cppclass SmartQuadtree[T]:
        cppclass const_iterator:
            const_iterator()
            T operator*()
            const_iterator operator++()
            bint operator==(const_iterator)
            bint operator!=(const_iterator)
            cpp_vector[const T*].const_iterator forward_begin()
            cpp_vector[const T*].const_iterator forward_end()
        cppclass iterator:
            iterator()
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        SmartQuadtree(double, double, double, double, unsigned int)
        cppbool insert(T)
        iterator begin()
        iterator end()
        const_iterator const_begin "begin" ()
        const_iterator const_end "end" ()

cdef extern from "cython_additions.hpp":
    double BoundaryXY_getX(long&)
    double BoundaryXY_getY(long&)
    double size_limit
    SmartQuadtree[long].iterator masked_begin(SmartQuadtree[long], PolygonMask*)
    SmartQuadtree[long].const_iterator masked_const_begin(SmartQuadtree[long], PolygonMask*)

cdef double BoundaryXY_getX(long& p):
    cdef object o = <object> (<void*> p)
    if isinstance(o, tuple) or isinstance(o, list):
        return <double> o[0]
    assert ("get_x" in dir(o))
    return <double> o.get_x()

cdef double BoundaryXY_getY(long& p):
    cdef object o = <object> (<void*> p)
    if isinstance(o, tuple) or isinstance(o, list):
        return <double> o[1]
    assert ("get_y" in dir(o))
    return <double> o.get_y()

cdef class Quadtree(object):
    """ Main class for quadtrees
    You must provide x-y coordinates for a center, x-y dimensions for width
    and height of the quadtree and a maximum size for data in each cell before
    subdivision (default: 8).

    For a quadtree centered on (0, 0) and ranging on [-10, 10]:
    >>> q = Quadtree(0., 0., 10., 10.)

    All methods attached to the Quadtree class are detailed below.
    """
    cdef PolygonMask* p
    cdef SmartQuadtree[long]* q

    def __cinit__(self, x0=None, y0=None, dim_x=None, dim_y=None, depth=8):
        self.p = NULL
        self.q = NULL
        if any([p is None for p in [x0, y0, dim_x, dim_y]]):
            print self.__doc__
            raise SyntaxError
        self.q = new SmartQuadtree[long](x0, y0, dim_x, dim_y, depth)

    def __dealloc__(self):
        if self.q != NULL:
            for x in self.elements():
                Py_DECREF(<object> (<void*> x))
            del self.q
        if self.p != NULL:
            del self.p

    def __repr__(self):
        cdef long count
        size_total = self.size()
        size_mask = self.size(False)
        string = "<smartquadtree.Quadtree at 0x%x>\n" % (<long> self.q)
        string += "Total number of elements: %ld\n" % size_total
        if (self.p != NULL):
            string += "Total number of elements inside mask: %ld\n"\
                % size_mask
            if (size_mask > 0):
                string += "First elements inside the mask:\n    "
        else:
            string += "No mask set\n"
            if (size_total > 0):
                string += "First elements:\n    "
        count = 0
        for i in self.elements():
            if (count > 3): return
            if (count == 3):
                string += "...  "
                break
            else: string += i.__repr__() + ",\n    "
            count += 1
        if (size_mask > 0 or size_total > 0):
            string = string[:-2]
        return string

    def set_mask(self, coords):
        """ Sets a polygon mask for future iterate functions.

        Just provide an iterable with (x, y) coordinates. Note that you do
        not need to close the polygon.
        >>> q.set_mask([ (0, 0), (1.5, 3), (3, 0) ])

        You can set the mask to None in order to reset it.
        >>> q.set_mask(None)

        """
        cdef vector[double] vec_x, vec_y
        if self.p != NULL:
            del self.p
            self.p = NULL
        if coords is None:
            return
        for (x, y) in coords:
            vec_x.push_back(x)
            vec_y.push_back(y)
        self.p = new PolygonMask(vec_x, vec_y, len(coords))

    def set_limitation(self, double size):
        """ Provides a criteria for stopping subdivisions.

        (see also: neighbours())

        The iterate_pair function provides a way to iterate over couples of
        items in a neighbourhood defined by the current cell and
        neighbouring cells. The set_limitation function lets you provide a
        minimum size for each cell of the quadtree.

        The quadtree will therefore stop subdivising cells when the size is
        reached.

        >>> q.set_limitation(2.)
        """
        global size_limit
        size_limit = size

    def insert(self, elt):
        """ Inserts an element into the quadtree.

        Just push in any element in the quadtree. The first element you push
        in determines what wil be pushed next.

        If the element is a tuple or a list, x coordinate is elt[0] et y
        coordinate is elt[1]. Otherwise the quadtree will try to access
        elt.get_x() and elt.get_y().

        It is not possible to switch between elt[0]/elt[1] and
        elt.get_x()/elt.get_y() after the first element has been inserted.

        All following examples are valid (provided that for q3, class Point
        implements get_x() and get_y().
        >>> q1.insert((1, 2))
        >>> q2.insert([1, 2])
        >>> q3.insert(Point(1, 2))

        """
        if self.q.insert(<long>(<void*> elt)):
            Py_INCREF(elt)

    def size(self, ignore_mask = True):
        """ Yields the size of the quadtree.

        This function is equivalent to iterate with a lambda incrementing a
        global variable.

        By default, any masks are ignored.
        >>> q.size()

        If a mask is set, you can provide False to stop ignoring masks.
        >>> q.size(False)

        """
        return sum(1 for x in self.elements(ignore_mask))

    def elements(self, ignore_mask = False):
        """ Iterates over elements in the quadtree (generator).

        You can iterate over all elements of the quadtree, and limit your
        parsing to elements inside a polygon defined in `set_mask`.

        As this method is a generator, elements are produced one by one and
        their position is adjusted in the quadtree after you ask for the next
        element. You are however certain never to get twice the same element.

        >>> for x in q.elements():
        >>>     print (x)

        If you want to print elements inside the following triangle:
        >>> q.set_mask([ (0, 0), (1.5, 3), (3, 0) ])
        >>> for x in q.elements():
        >>>     print (x)

        You can ignore the mask with the `ignore_mask` parameter
        (default: False)
        >>> for x in q.elements(ignore_mask = True):
        >>>     print (x)
        """
        cdef SmartQuadtree[long].iterator it
        if self.p is NULL or ignore_mask:
            it = self.q.begin()
        else:
            it = masked_begin(deref(self.q), self.p)
        while (it != self.q.end()):
            yield (<object>(<void*> deref(it)))
            inc(it)

    def neighbour_elements(self, ignore_mask = False):
        """ Iterates over pair of neighbours in the quadtree (generator).

        You can iterate over all pairs of elements of the quadtree, and limit
        your parsing to elements inside a polygon defined in `set_mask`.

        As this method is a generator, elements are produced one by one. Note
        that if (a, b) is produced by the generator, (b, a) will not be yielded.

        Be warned that if elements are moved through this iterator, the quadtree
        will be obsolete. Therefore it is considered unsafe to use this iterator
        to move the positions of elements.


        >>> for a, b in q.neighbour_elements():
        >>>     print (a, b)

        If you want to print elements inside the following triangle:
        >>> q.set_mask([ (0, 0), (1.5, 3), (3, 0) ])
        >>> for a, b in q.neighbour_elements():
        >>>     print (a, b)

        You can ignore the mask with the `ignore_mask` parameter
        (default: False)
        >>> for a, b in q.neighbour_elements(ignore_mask = True):
        >>>     print (a, b)
        """
        cdef SmartQuadtree[long].const_iterator it
        if self.p is NULL or ignore_mask:
            it = self.q.const_begin()
        else:
            it = masked_const_begin(deref(self.q), self.p)
        cdef cpp_vector[const long*].const_iterator j
        while (it != self.q.const_end()):
            j = it.forward_begin()
            while (j != it.forward_end()):
                yield(<object>(<void*> deref(it)),
                        <object>(<void*> deref(deref(j))))
                inc(j)
            inc(it)
