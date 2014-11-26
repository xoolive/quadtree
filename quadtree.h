#ifndef QUADTREE_H
#define QUADTREE_H

#include <cassert>
#include <cstdlib>

#include <list>
#include <vector>
#include <iostream>
#include <unordered_map>

#include "neighbour.h"

class ExtendedQuadtree;
class Boundary;

class PolygonMask
{
private:
  //! Nb of vertices in the polygon
  int size;

  //! Coordinates of the vertices of the polygon
  std::vector<float> polyX, polyY;

  //! Auxiliary variables
  std::vector<float> constant, multiple;

  // see http://alienryderflex.com/polygon/
  void precompute();

public:

  //! Constructor
  PolygonMask(std::vector<float> x, std::vector<float> y, int size);

  //! Return the number of vertices of the polygon
  int getSize() const { return size; }

  //! Returns whether a point of coordinates (x, y) is inside the polygon
  // see http://alienryderflex.com/polygon/
  bool pointInPolygon(float x, float y) const;

  //! Returns a different polygon mask clipped by the boundary box
  // see http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
  PolygonMask clip(const Boundary& box) const;

  friend class Test_PolygonMask;

};

class Boundary
{
  //! Coordinates for the center of the box
  float center_x, center_y ;

  //! Dimension from the center to the border of the box
  float dim_x, dim_y;

  //! Find the coordinates of the data
  float (*x)(void*), (*y)(void*);

  //! Give some limitation to the size of the cells. This function shall return
  //! true if the the Boundary box must not be subdivided.
  bool (*limitfct)(Boundary*);

  //! Store the result of limitation in order to avoid recomputation
  bool limit;

  //! Limitation function; returns false if unset
  bool limitation() {
    if (limitfct == NULL) return false;
    bool limit = limitfct(this);
    return limit;
  }

  //! Return the number of points covered by the polygon mask m
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

  //! Returns the x-coordinate of the center of the boundary box
  inline float getX() const { return center_x; }

  //! Returns the y-coordinate of the center of the boundary box
  inline float getY() const { return center_y; }

  //! Returns the x dimension of the boundary box
  inline float getDimX() const { return dim_x; }

  //! Returns the y dimension of the boundary box
  inline float getDimY() const { return dim_y; }

  //! L1 norm
  inline float norm_l1() const { return (dim_x + dim_y); }

  //! Max distance between two data in the box
  inline float norm_infty() const { return (dim_x < dim_y ? dim_x : dim_y); }

  //! Returns true if the coordinates are left of the boundary box
  bool leftOf(float x, float y) const { return (x < center_x-dim_x-1e-4); }

  //! Returns true if the coordinates are right of the boundary box
  bool rightOf(float x, float y) const { return (x > center_x+dim_x+1e-4); }

  //! Returns true if the coordinates are bottom of the boundary box
  bool bottomOf(float x, float y) const { return (y < center_y-dim_y-1e-4); }

  //! Returns true if the coordinates are up of the boundary box
  bool upOf(float x, float y) const { return (y > center_y+dim_y+1e-4); }

  //! Returns the intersection with the left boundary of the box
  void interLeft(float, float, float, float, float&, float&);

  //! Returns the intersection with the right boundary of the box
  void interRight(float, float, float, float, float&, float&);

  //! Returns the intersection with the bottom boundary of the box
  void interBottom(float, float, float, float, float&, float&);

  //! Returns the intersection with the up boundary of the box
  void interUp(float, float, float, float, float&, float&);

  //! Types of methods leftOf, rightOf, bottomOf, upOf
  typedef bool (Boundary::*OUTSIDE_TEST) (float, float) const;

  //! Types of methods interLeft, interRight, interBottom, interUp
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

  // Just to help not iterating twice
  std::list<void*> already;

  // We keep a map of who is where
  std::unordered_map<void*, ExtendedQuadtree*> where;

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

  //! Destructor
  ~ExtendedQuadtree();

  //! Find same level neighbour in determined direction
  ExtendedQuadtree* samelevel(unsigned char) const;

  //! Insert one piece of data to the quadrant
  //! Returns true if the data has been inserted
  bool insert(void*);

  //! Update a data in current subtree, returns true if changed cell
  bool updateData(void* p);

  //! Removes a data in current subtree
  void removeData(void* p);

  //! Returns true if the current cell may contain the data
  bool contains(void* p) { return b.contains(p); }

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

  //! Iterate something for all items
  void iterate(bool (*apply)(void*));

  //! Iterate something for all pairs of items
  void iterate(void (*apply)(void*, void*));

  //! Iterate some function for all pairs of items, vectorised version
  void iterateby4(void (*apply)(void*, void*), void (*applyby4)(void*,void**));

  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);
  friend class Test_ExtendedQuadtree;

};

#endif // QUADTREE_H
