#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <vector>
#include <iostream>

enum Direction {
  EAST, NORTHEAST, NORTH, NORTHWEST, WEST, SOUTHWEST, SOUTH, SOUTHEAST
};


class Neighbour
{

  //! Max level generated for _x and _y
  int maxlevel;

  //! Basic movement backwards in x
  unsigned int _x;

  //! Basic movement backwards in y
  unsigned int _y;

  unsigned int east, northeast, north, northwest,
               west, southwest, south, southeast;

  Neighbour() { maxlevel = 0; _x = 0; _y = 0; }

  static Neighbour* getInstance() {
    static Neighbour* instance = NULL;
    if (NULL == instance) {
      instance            = new Neighbour();
      instance->east      = 1;
      instance->north     = 2;
      instance->northeast = 3;
    }
    return instance;
  }

  static void generate(Neighbour* instance, unsigned int level)
  {
    for (int i = instance->maxlevel; i<=level; ++i)
    {
      instance->_x = (instance->_x << 2) + 1;
      instance->_y = (instance->_y << 2) + 2;
    }
    instance->maxlevel  = level;
    instance->west      = instance->_x;
    instance->south     = instance->_y;
    instance->northwest = instance->north + instance->west;
    instance->southwest = instance->south + instance->west;
    instance->southeast = instance->south + instance->east;
  }

  static unsigned int direction(Neighbour* instance, enum Direction dir)
  {
    switch (dir) {
    case EAST: return instance->east;
    case NORTHEAST: return instance->northeast;
    case NORTH: return instance->north;
    case NORTHWEST: return instance->northwest;
    case WEST: return instance->west;
    case SOUTHWEST: return instance->southwest;
    case SOUTH: return instance->south;
    default: return instance->southeast;
    }
  }

public:

  static unsigned int samelevel(unsigned int x, enum Direction dir,
                                unsigned int level)
  {
    Neighbour* instance = getInstance();
    if (level > instance->maxlevel) generate(instance, level);
    unsigned int tx = instance->_x/*[instance->maxlevel+1]*/;
    unsigned int ty = instance->_y/*[instance->maxlevel+1]*/;
    unsigned int d = direction(instance, dir);
    return (((x | ty) + (d & tx)) & tx) | (((x | tx) + (d & ty)) & ty);
  }


};

#endif // NEIGHBOUR_H
