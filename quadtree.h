#ifndef QUADTREE_H
#define QUADTREE_H

#include <list>
#include <cassert>
#include <cstdlib> // NULL
#include <iostream>
#include <vector>

#include "neighbour.h"

class ExtendedQuadtree;
class Boundary;

class PolygonMask
{
private:
  int size;
  std::vector<float> polyX, polyY;
  std::vector<float> constant, multiple;

  // see http://alienryderflex.com/polygon/
  void precompute();

public:

  PolygonMask(std::vector<float> x, std::vector<float> y, int size);

  int getSize() const { return size; }

  // see http://alienryderflex.com/polygon/
  bool pointInPolygon(float x, float y) const;

  // see http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
  PolygonMask clip(const Boundary& box) const;

  friend class Test_PolygonMask;

};

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

  int coveredByPolygon(const PolygonMask& m);

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

  bool leftOf(float x, float y) const { return (x < center_x - dim_x - 1e-4); }
  bool rightOf(float x, float y) const { return (x > center_x + dim_x + 1e-4); }
  bool bottomOf(float x, float y) const { return (y < center_y - dim_y - 1e-4); }
  bool upOf(float x, float y) const { return (y > center_y + dim_y + 1e-4); }

  void interLeft(float, float, float, float, float&, float&);
  void interRight(float, float, float, float, float&, float&);
  void interBottom(float, float, float, float, float&, float&);
  void interUp(float, float, float, float, float&, float&);

  typedef bool (Boundary::*OUTSIDE_TEST) (float, float) const;
  typedef bool (Boundary::*INTERSECT)
    (float, float, float, float, float&, float&) const;

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

  // Directions corresponding to children nodes
  static const unsigned char diags[4];

  // Data attached to the quadrant
  std::list<void*> points;

  // Capacity of each cell
  const unsigned int capacity;

  // Ancestor
  ExtendedQuadtree* ancestor;

  //! Increments the delta in direction dir
  //! Returns true if you have children
  bool incrementDelta(unsigned char dir, bool flag = true);

  //! Update the no more non reflexive delta
  void updateDiagonal(unsigned char diagdir, unsigned char dir, int delta);

  //! Updates the delta in direction dir if neighbour of same level has
  //! children nodes
  void updateDelta(unsigned char dir);

public:

  //! Constructor
  ExtendedQuadtree(float center_x, float center_y, float dim_x, float dim_y,
                   unsigned int capacity) :
    b(center_x, center_y, dim_x, dim_y), location(0), level(0),
    capacity(capacity), ancestor(this)
  {
    children[0] = NULL; children[1] = NULL;
    children[2] = NULL; children[3] = NULL;
    delta[0] = 2; delta[1] = 2; delta[2] = 2; delta[3] = 2;
    delta[4] = 2; delta[5] = 2; delta[6] = 2; delta[7] = 2;
  }

  //! Constructor of a child quadtree
  // SW -> 0, SE -> 1, NW -> 2, NE -> 3
  ExtendedQuadtree(const ExtendedQuadtree&, unsigned char);

  //! Find same level neighbour in determined direction
  ExtendedQuadtree* samelevel(unsigned char) const;

  //! Insert one piece of data to the quadrant
  bool insert(void*);

  //! Returns the subquadrant pointed by location code
  ExtendedQuadtree* getQuadrant(unsigned long location,
                                unsigned short level) const;

  //! Returns the data embedded to current quadrant
  inline const std::list<void*>& getPoints() const { return points; }

  //! Returns a point to the proper child 0->SW, 1->SE, 2->NW, 3->NE
  inline const ExtendedQuadtree* getChild(unsigned char i) const
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
  inline unsigned long getLocation() const { return location; }

  //! Get the level, only for debugging purposes
  inline unsigned char getLevel() const { return level; }

  //! Get the max size of children data lists
  unsigned long getDataSize() const;

  //! Get the depth of the quadtree
  unsigned char getDepth() const;

  //! Iterate something for all items
  //! Adjusts the quadtree if the items move
  void iterate(const PolygonMask& m, bool (*apply)(void*));
  void iterate(bool (*apply)(void*));

  //! Iterate something for all pairs of items
  void iterate(void (*apply)(void*, void*));

  //! Iterate some function for all pairs of items, vectorised version
  void iterateby4(void (*apply)(void*, void*), void (*applyby4)(void*,void**));

  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);
  friend class Test_ExtendedQuadtree;

};

#endif // QUADTREE_H
