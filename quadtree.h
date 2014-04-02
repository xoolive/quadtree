#ifndef QUADTREE_H
#define QUADTREE_H

#include <list>
#include <cassert>
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

  // Store the result of limitation in order to avoid recomputation
  bool limit;

  bool limitation() {
    if (limitfct == NULL) return false;
    bool limit = limitfct(this);
    return limit;
  }

public:

  //! Default constructor
  Boundary(float cx, float cy, float dx, float dy) :
    center_x(cx), center_y(cy), dim_x(dx), dim_y(dy),
    x(NULL), y(NULL), limitfct(NULL) {};

  //! Is this point included in the box?
  bool contains(float x, float y);

  //! Is this data included in the box?
  inline bool contains(void* pt) { return contains(x(pt), y(pt)); }

  //! L1 norm
  inline float norm_l1() { return (dim_x + dim_y); }

  //! Max distance between two data in the box
  inline float norm_infty() { return (dim_x < dim_y ? dim_x : dim_y); }

  friend class ExtendedQuadtree;
  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);

};

class Test_ExtendedQuadtree;
class ExtendedQuadtree
{
  // Delimitates the quadrant
  Boundary b;

  // Binary representation of the location code
  std::size_t location;

  // Current level of the quadrant
  int level;

  // Level differences with the neighbours
  // 0: adjacent quadrant is of same level
  // 1: adjacent quadrant is of larger level (i.e. smaller in size)
  // 2: adjacent quadrant is out of area
  // 3: adjacency is not reflexive (diagonal)
  // -1, ...: adjacent quadrant is of smaller level by -n (i.e. larger in size)
  int delta[8];

  // Children nodes
  ExtendedQuadtree *children[4];

  // Data attached to the quadrant
  std::list<void*> points;

  // Capacity of each cell
  const int capacity;

  // Ancestor
  ExtendedQuadtree* ancestor;

  //! Increments the delta in direction dir
  //! Returns true if you have children
  bool incrementDelta(unsigned int dir, bool flag = true);

  //! Update the no more non reflexive delta
  void updateDiagonal(unsigned int diagdir, unsigned int dir, int delta);

  //! Updates the delta in direction dir if neighbour of same level has
  //! children nodes
  void updateDelta(unsigned int dir);

public:

  //! Constructor
  ExtendedQuadtree(float center_x, float center_y, float dim_x, float dim_y,
                   int capacity) :
    b(center_x, center_y, dim_x, dim_y), location(0), level(0),
    capacity(capacity), ancestor(this)
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
  inline const std::list<void*>& getPoints() const { return points; }

  //! Returns a point to the proper child 0->SW, 1->SE, 2->NW, 3->NE
  inline const ExtendedQuadtree* getChild(unsigned int i) const
  {
    assert (i<4);
    return children[i];
  }

  //! Sets the getters in the Boundary
  inline void setXYFcts(float (*x)(void*), float (*y)(void*))
  { b.x = x; b.y = y; }

  //! Sets the limitation function; if NULL, behaves as no restriction
  inline void setLimitation(bool (*limit)(Boundary*))
  {
    assert(limit != NULL);
    b.limitfct = limit;
    b.limit = b.limitation();
  }

  //! Get the location
  inline std::size_t getLocation() const { return location; }

  //! Get the level, only for debugging purposes
  inline std::size_t getLevel() const { return level; }

  //! Get the max size of children data lists
  std::size_t getDataSize() const;

  //! Get the depth of the quadtree
  std::size_t getDepth() const;

  //! Iterate something for all items
  //! Adjusts the quadtree if the items move
  void iterate(bool (*apply)(void*));

  //! Iterate something for all pairs of items
  void iterate(void (*apply)(void*, void*));

  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);
  friend class Test_ExtendedQuadtree;

};

#endif // QUADTREE_H
