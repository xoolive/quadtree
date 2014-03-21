#ifndef QUADTREE_H
#define QUADTREE_H

#include <list>
#include <cstddef> // size_t
#include <cstdlib> // NULL
#include <iostream>

#include "neighbour.h"

class ExtendedQuadtree;

class Boundary
{
  // Coordinates for the center of the box
  float center_x;
  float center_y;

  // Dimension from the center to the border of the box
  float dim_x;
  float dim_y;

  // Find the coordinates of the data
  float (*x)(void*);
  float (*y)(void*);

  // Give some limitation to the size of the cells
  bool (*limitfct)(Boundary*);

  bool limitation() {
    if (limitfct == NULL) return false;
    return limitfct(this);
  }

public:

  //! Default constructor
  Boundary(float cx, float cy, float dx, float dy) :
    center_x(cx), center_y(cy), dim_x(dx), dim_y(dy),
    x(NULL), y(NULL), limitfct(NULL) {};

  //! Is this point included in the box?
  bool contains(float x, float y);

  //! Is this data included in the box?
  bool contains(void* pt) { return contains(x(pt), y(pt)); }

  //! Max distance between two data in the box
  float norm_l1() { return (dim_x + dim_y); }

  //! Max distance between two data in the box
  float sq_norm_l2() { return (dim_x * dim_x + dim_y * dim_y); }

  friend class ExtendedQuadtree;
  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);

};

class ExtendedQuadtree
{
  // Delimitates the quadrant
  Boundary b;

  // Binary representation of the location code
  std::size_t location;

  // Current level of the quadrant
  int level;

public: // TODO remove public
  // Level differences with the neighbours
  // 0: adjacent quadrant is of same level
  // 1: adjacent quadrant is of larger level (i.e. smaller in size)
  // 2: adjacent quadrant is out of area
  // -1, ...: adjacent quadrant is of smaller level by -n (i.e. larger in size)
  int delta[8];

  // Children nodes
  ExtendedQuadtree *children[4];

private:
  // Data attached to the quadrant
  std::list<void*> points;

  // Capacity of each cell
  const int capacity;

  // Ancestor
  ExtendedQuadtree* ancestor;

  //! Increments the delta in direction dir
  //! Returns true if you have children
  bool incrementDelta(unsigned int dir);

  //! Updates the delta in direction dir if neighbour of same level has
  //! children nodes
  void updateDelta(unsigned int dir);

public:

  //! Constructor
  ExtendedQuadtree(float center_x, float center_y, float dim_x, float dim_y,
                   int capacity) :
    b(center_x, center_y, dim_x, dim_y), location(0), level(0),
    /*dn(2), dw(2), ds(2), de(2),*/ capacity(capacity), ancestor(this)
  {
    for (int i = 0; i<4; ++i) children[i] = NULL;
    for (int i = 0; i<8; ++i) delta[i] = 2;
  }

  //! Constructor of a child quadtree
  // SW -> 0, SE -> 1, NW -> 2, NE -> 3
  ExtendedQuadtree(const ExtendedQuadtree&, int);

  //! Find same level neighbour in determined direction
  ExtendedQuadtree* samelevel(unsigned int) const;

  //! Insert one piece of data to the quadrant
  bool insert(void*);

  //! Returns the subquadrant pointed by location code
  ExtendedQuadtree* getQuadrant(std::size_t location, std::size_t level) const;

  //! Returns the data embedded to current quadrant
  const std::list<void*>& getPoints() const { return points; }

  //! Sets the getters in the Boundary
  void setXYFcts(float (*x)(void*), float (*y)(void*)) { b.x = x; b.y = y; }

  //! Sets the limitation function; if NULL, behaves as no restriction
  void setLimitation(bool (*limit)(Boundary*)) { b.limitfct = limit; }

  //! Get the location
  std::size_t getLocation() const { return location; }

  //! Get the level, only for debugging purposes
  std::size_t getLevel() const { return level; }

  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);

};

#endif // QUADTREE_H
