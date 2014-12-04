# distutils: language = c++
# -*- coding: utf-8 -*-

from libcpp cimport bool as cppbool
from libcpp.vector cimport vector
from cpython.ref cimport Py_INCREF, Py_DECREF

cdef extern from "quadtree.h":
    cdef cppclass PolygonMask:
        PolygonMask(vector[float], vector[float], int)
    cdef cppclass SmartQuadtree:
        SmartQuadtree(float, float, float, float, unsigned int)
        void setXYFcts(float (*x)(void*), float (*y)(void*))
        cppbool insert(void*)
        void iterate(cppbool (*fct)(void*))
        void iterate(PolygonMask, cppbool (*fct)(void*))

cdef object global_lbd

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

cdef cppbool apply_lambda(void* p):
    global global_lbd
    return <bint> global_lbd(<object> (p))

cdef cppbool decref_all(void* p):
    Py_DECREF(<object> p)
    return <bint> False

cdef class Quadtree(object):

    cdef PolygonMask* p
    cdef SmartQuadtree* q
    cdef object is_tuple_or_list

    def __cinit__(self, x0, y0, dim_x, dim_y, depth):
        self.p = NULL
        self.q = new SmartQuadtree(x0, y0, dim_x, dim_y, depth)
        self.is_tuple_or_list = None

    def __dealloc__(self):
        self.q.iterate(decref_all)
        del self.q
        if self.p != NULL:
            del self.p

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
        if coords is None:
            return
        for (x, y) in coords:
            vec_x.push_back(x)
            vec_y.push_back(y)
        self.p = new PolygonMask(vec_x, vec_y, len(coords))

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
        if self.p == NULL:
            self.q.iterate(apply_lambda)
        else:
            self.q.iterate(self.p[0], apply_lambda)




