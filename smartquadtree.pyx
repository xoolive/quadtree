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

cdef extern from "quadtree.h":
    cdef cppclass Boundary:
        float norm_infty()
    cdef cppclass PolygonMask:
        PolygonMask(vector[float], vector[float], int)
    cdef cppclass SmartQuadtree:
        SmartQuadtree(float, float, float, float, unsigned int)
        void setXYFcts(float (*x)(void*), float (*y)(void*))
        void setLimitation(cppbool (*limit)(Boundary*))
        cppbool insert(void*)
        void iterate(cppbool (*fct)(void*))
        void iterate(void (*fct)(void*, void*))
        void iterate(PolygonMask, cppbool (*fct)(void*))

cdef object global_lbd, global_flag
cdef float size_limit
cdef object string

cdef float c_getx(void* p):
    cdef object o = <object> (p)
    assert ("get_x" in dir(o))
    return <float> o.get_x()

cdef float c_getidx0(void* p):
    cdef object o = <object> (p)
    assert (isinstance(o, tuple) or isinstance(o, list))
    return <float> o[0]

cdef float c_gety(void* p):
    cdef object o = <object> (p)
    assert ("get_y" in dir(o))
    return <float> o.get_y()

cdef float c_getidx1(void* p):
    cdef object o = <object> (p)
    assert (isinstance(o, tuple) or isinstance(o, list))
    return <float> o[1]


class movable_elt(object):
    """ Decorator for callbacks that modify x-y coordinates of elements. """

    def __init__(self, f):
        self.func = f

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def __repr__(self):
        return self.func.__doc__


class static_elt(object):
    """ Decorator for callbacks that do not modify x-y coordinates of elements. """

    def __init__(self, f):
        self.func = f

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def __repr__(self):
        return self.func.__doc__


cdef cppbool apply_lambda(void* p):
    global global_lbd, global_flag
    global_lbd(<object> (p))
    return (<bint> global_flag)

cdef void apply_lambda_pair(void* p1, void* p2):
    global global_lbd
    global_lbd(<object> (p1), <object> (p2))
    return

cdef cppbool decref_all(void* p):
    Py_DECREF(<object> p)
    return <bint> False

cdef cppbool limitation(Boundary* b):
    global size_limit
    cdef float sq_size = b.norm_infty()
    return <bint> (sq_size < size_limit)

cdef class Quadtree(object):
    """ Main class for quadtrees
    You must provide x-y coordinates for a center, x-y dimensions for width
    and height of the quadtree and a maximum depth for subdivisions.

    For a quadtree centered on (0, 0) and ranging on [-10, 10]:
    >>> q = Quadtree(0., 0., 10., 10., 16)

    All methods attached to the Quadtree class are detailed below.
    """
    cdef PolygonMask* p
    cdef SmartQuadtree* q
    cdef object is_tuple_or_list
    cdef long count

    def __cinit__(self, x0, y0, dim_x, dim_y, depth):
        self.p = NULL
        self.q = new SmartQuadtree(x0, y0, dim_x, dim_y, depth)
        self.is_tuple_or_list = None
        self.count = 0

    def __dealloc__(self):
        self.q.iterate(decref_all)
        del self.q
        if self.p != NULL:
            del self.p

    def __repr__(self):
        global string
        size_total = self.size()
        size_mask = self.size(False)
        string = "<smartquadtree.Quadtree at 0x%x>\n" % (<long> self.q)
        string += "Total number of elements: %ld\n" % size_total
        if (self.p != NULL):
            string += "Total number of elements inside mask: %ld\n"\
                % size_mask
            if (size_mask > 0):
                string += "First elements inside the mask:\n"
        else:
            string += "No mask set\n"
            if (size_total > 0):
                string += "First elements:\n    "
        self.count = 0
        @static_elt
        def print_first(p):
            global string
            if (self.count > 3): return
            if (self.count == 3):
                string += "...  "
                self.count += 1
                return
            string += p.__repr__() + ", "
            self.count += 1
        self.iterate(print_first)
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
        cdef vector[float] vec_x, vec_y
        if self.p != NULL:
            del self.p
            self.p = NULL
        if coords is None:
            return
        for (x, y) in coords:
            vec_x.push_back(x)
            vec_y.push_back(y)
        self.p = new PolygonMask(vec_x, vec_y, len(coords))

    def set_limitation(self, size):
        """ Provides a criteria for stopping subdivisions.

        (see also: iterate_pair)

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
        self.q.setLimitation(limitation)

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
        if self.is_tuple_or_list is None:
            self.is_tuple_or_list =\
                isinstance(elt, tuple) or isinstance(elt, list)
            if self.is_tuple_or_list:
                self.q.setXYFcts(c_getidx0, c_getidx1)
            else:
                self.q.setXYFcts(c_getx, c_gety)
        Py_INCREF(elt)
        self.q.insert(<void*> elt)

    def iterate(self, lbd):
        """ Iterates over elements in the quadtree.

        You can provide a function or a lambda to be applied on all elements
        of the quadtree. If you provided a polygon beforehand (see
        set_mask), the function (resp. lambda) is only applied to elements
        within the polygon.

        If the function changes x-y coordinates of the element, it shall
        return True. If it does not change x-y coordinates, it shall return
        False.

        If the function does not return a Python boolean, it may be
        interpreted as False and the element may be incorrectly positioned
        inside the structure. (Better think undefined behaviour)

        >>> def f_print(p): print (p)
        >>> q.iterate(f_print)

        Also with lambdas:
        >>> q.iterate(lambda p: print (p))

        If you want to print elements inside the following triangle:
        >>> q.set_mask([ (0, 0), (1.5, 3), (3, 0) ])
        >>> q.iterate(lambda p: print (p))
        """
        global global_lbd
        global_lbd = lbd
        movable_flag = isinstance(global_lbd, movable_elt)
        static_flag = isinstance(global_lbd, static_elt)
        assert (movable_flag or static_flag),\
            "You must decorate your function with @static_elt or @movable_elt"
        global_flag = movable_flag
        if self.p == NULL:
            self.q.iterate(apply_lambda)
        else:
            self.q.iterate(self.p[0], apply_lambda)

    def iterate_pair(self, lbd):
        """ Iterates over pairs of neighbouring elements in the quadtree.

        (see also: set_limitation)

        The iterate_pair function provides a way to iterate over couples of
        items in a neighbourhood defined by the current cell and
        neighbouring cells. The set_limitation function lets you provide a
        minimum size for each cell of the quadtree, in order to NOT miss any
        neighbouring pair you would want to consider here.

        >>> def f_close(p1, p2): print ([p1, p2])
        >>> q.iterate_pair(f_close)

        """
        global global_lbd
        global_lbd = lbd
        self.q.iterate(apply_lambda_pair)

    def size(self, ignore_mask = True):
        """ Yields the size of the quadtree.

        This function is equivalent to iterate with a lambda incrementing a
        global variable.

        By default, any masks are ignored.
        >>> q.size()

        If a mask is set, you can provide False to stop ignoring masks.
        >>> q.size(False)

        """
        global global_lbd
        self.count = 0
        @static_elt
        def fun_count(p):
            self.count += 1
        global_lbd = fun_count
        if (ignore_mask or self.p == NULL):
            self.q.iterate(apply_lambda)
        else:
            self.q.iterate(self.p[0], apply_lambda)
        return self.count

