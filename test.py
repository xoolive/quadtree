from random import random
import smartquadtree


# Return nothing SEEMS TO BE like returning False
def print_cb(p):
    print (p)


class Point(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "(%f, %f)" % (self.x, self.y)

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

# Quadtree of objects
q1 = smartquadtree.Quadtree(0, 0, 10, 10, 16)
for a in range(20):
    q1.insert(Point(random()*20-10, random()*20-10))
q1.iterate(print_cb)  # works with functions

# Quadtree of array
q2 = smartquadtree.Quadtree(0, 0, 10, 10, 16)
for a in range(20):
    q2.insert([random()*20-10, random()*20-10])
q2.iterate(lambda p: print_cb(p))  # works with lambdas as well

print ("set mask")
q2.set_mask([(-5, -5), (-5, 5), (5, 5), (5, -5)])
q2.iterate(print_cb)
