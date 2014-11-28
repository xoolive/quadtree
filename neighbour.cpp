/*
 * Interface for a smart version of quadtrees specialised for tracking moving
 * objects. When iterating over elements in the quadtree, a simple flag
 * (boolean) indicates whether the element might have moved to a neighbouring
 * subdivision.
 *
 * Xavier Olive, 28 nov. 2014
 */

#include "neighbour.h"

int Neighbour::maxlevel = 0;
unsigned int Neighbour::_x = 0;
unsigned int Neighbour::_y = 0;
unsigned int Neighbour::directions[8] = {1, 3, 2, 0, 0, 0, 0, 0};

