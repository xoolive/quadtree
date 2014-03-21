#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

enum Direction {
  EAST, NORTHEAST, NORTH, NORTHWEST, WEST, SOUTHWEST, SOUTH, SOUTHEAST
};


class Neighbour
{
  //! Max level generated for _x and _y
  static int maxlevel;

  //! Basic movement backwards in x
  static unsigned int _x;

  //! Basic movement backwards in y
  static unsigned int _y;

  //! Codes for all direction movements
  static unsigned int directions[8];

  //! Generates codes for all direction movements
  static void generate(unsigned int level)
  {
    for (int i = maxlevel; i<=level; ++i)
    {
      _x = (_x << 2) + 1;
      _y = (_y << 2) + 2;
    }
    maxlevel  = level;
    directions[WEST]      = _x;
    directions[SOUTH]     = _y;
    directions[NORTHWEST] = directions[NORTH] + directions[WEST];
    directions[SOUTHWEST] = directions[SOUTH] + directions[WEST];
    directions[SOUTHEAST] = directions[SOUTH] + directions[EAST];
  }

public:

  //! Yields the location code for the neighbour of same level in direction dir
  static unsigned int samelevel(unsigned int x, unsigned int dir,
                                unsigned int level)
  {
    if (level > maxlevel) generate(level);
    return (((x | _y) + (directions[dir] & _x)) & _x) |
      (((x | _x) + (directions[dir] & _y)) & _y);
  }
};

#endif // NEIGHBOUR_H
