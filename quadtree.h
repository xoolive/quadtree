/*
 * Interface for a smart version of quadtrees specialised for tracking moving
 * objects. Iterators and their const counterparts let you process objects and
 * ajust the quadtree structure when points are moved.
 *
 * Xavier Olive, 28 nov. 2014
 */

#ifndef QUADTREE_H
#define QUADTREE_H

#include <cassert>
#include <cstdlib>

#include <list>
#include <vector>
#include <iostream>
#include <unordered_map>

#include "neighbour.h"

template<class T> class SmartQuadtree;
template<class T> class MaskedQuadtree;

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

template<typename T>
struct BoundaryXY
{
  static double getX(const T& p) { return p.x; }
  static double getY(const T& p) { return p.y; }
};

template<typename T>
struct BoundaryLimit
{
  static bool limitation(const T& box) { return false; }
};

class Boundary
{
  //! Coordinates for the center of the box
  float center_x, center_y ;

  //! Dimension from the center to the border of the box
  float dim_x, dim_y;

  //! Store the result of limitation in order to avoid recomputation
  bool limit;

  //! Return the number of points covered by the polygon mask m
  int coveredByPolygon(const PolygonMask& m) const;

public:

  //! Default constructor
  Boundary(float cx, float cy, float dx, float dy) :
    center_x(cx), center_y(cy), dim_x(dx), dim_y(dy), limit(false) {};

  Boundary(const Boundary& b) :
    center_x(b.center_x), center_y(b.center_y),
    dim_x(b.dim_x), dim_y(b.dim_y), limit(false) {};

  //! Is this point included in the box?
  bool contains(float x, float y);

  //! Is this data included in the box?
  template<typename T>
  inline bool contains(const T& pt)
  { return contains(BoundaryXY<T>::getX(pt), BoundaryXY<T>::getY(pt)); }

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

  template<typename T> friend class SmartQuadtree;
  template<typename T>
  friend std::ostream& operator<< (std::ostream&, const SmartQuadtree<T>&);

};

class Test_SmartQuadtree;

// Necessary on some compilers in order to be befriended
template<typename T>
std::ostream& operator<< (std::ostream&, const SmartQuadtree<T>&);

template<typename T>
class SmartQuadtree
{
  // Delimitates the quadrant
  Boundary b;

  // Binary representation of the location code
  std::size_t location;

  // Current level of the quadrant
  std::size_t level;

  // Level differences with the neighbours
  // 0: adjacent quadrant is of same level
  // 1: adjacent quadrant is of larger level (i.e. smaller in size)
  // 2: adjacent quadrant is out of area
  // 3: adjacency is not reflexive (diagonal)
  // -1, ...: adjacent quadrant is of smaller level by -n (i.e. larger in size)
  int delta[8];

  // Children nodes
  SmartQuadtree<T> *children[4];

  // Directions corresponding to children nodes
  static const unsigned char diags[4];

  // Data attached to the quadrant
  std::list<T> points;

  // We keep a map of who is where
  std::unordered_map<T*, SmartQuadtree*> where;

  // All leaves of the Quadtree, in order
  std::list<SmartQuadtree*> leaves;

  // Capacity of each cell
  const unsigned int capacity;

  // Ancestor
  SmartQuadtree<T>* ancestor;

  //! Increments the delta in direction dir
  //! Returns true if you have children
  bool incrementDelta(unsigned char dir, bool flag = true);

  //! Update the no more non reflexive delta
  void updateDiagonal(unsigned char diagdir, unsigned char dir, int delta);

  //! Updates the delta in direction dir if neighbour of same level has
  //! children nodes
  void updateDelta(unsigned char dir);

public:

  struct const_iterator;
  struct iterator;

  //! Constructor
  SmartQuadtree<T>(float center_x, float center_y, float dim_x, float dim_y,
                unsigned int capacity) :
    b(center_x, center_y, dim_x, dim_y), location(0), level(0),
    capacity(capacity), ancestor(this)
  {
    children[0] = NULL; children[1] = NULL;
    children[2] = NULL; children[3] = NULL;
    delta[0] = 2; delta[1] = 2; delta[2] = 2; delta[3] = 2;
    delta[4] = 2; delta[5] = 2; delta[6] = 2; delta[7] = 2;
    leaves.push_back(this);
  }

  //! Constructor of a child quadtree
  // SW -> 0, SE -> 1, NW -> 2, NE -> 3
  SmartQuadtree<T>(const SmartQuadtree<T>&, unsigned char,
                   typename std::list<SmartQuadtree<T>*>::iterator&);

  //! Destructor
  ~SmartQuadtree<T>();

  //! Iterator
  SmartQuadtree<T>::iterator begin();

  //! Iterator
  SmartQuadtree<T>::iterator end();

  //! Iterator (const version)
  SmartQuadtree<T>::const_iterator begin() const;

  //! Iterator (const version)
  SmartQuadtree<T>::const_iterator end() const;

  //! Find same level neighbour in determined direction
  SmartQuadtree<T>* samelevel(unsigned char) const;

  //! Insert one piece of data to the quadrant
  //! Returns true if the data has been inserted
  bool insert(T);

  //! Update a data in current subtree, returns true if changed cell
//   bool updateData(T& p);

  //! Removes a data in current subtree
  void removeData(T& p);

  //! Returns true if the current cell may contain the data
  bool contains(T& p) { return b.contains(&p); }

  //! Returns the subquadrant pointed by location code
  SmartQuadtree<T>* getQuadrant(unsigned long location,
                                unsigned short level) const;

  //! Returns the data embedded to current quadrant
  inline const std::list<T>& getPoints() const { return points; }

  //! Returns a point to the proper child 0->SW, 1->SE, 2->NW, 3->NE
  inline const SmartQuadtree<T>* getChild(unsigned char i) const
  {
    assert (i<4);
    return children[i];
  }

  //! Get the location
  inline unsigned long getLocation() const { return location; }

  //! Get the level, only for debugging purposes
  inline unsigned char getLevel() const { return level; }

  //! Get the max size of children data lists
  unsigned long getDataSize() const;

  //! Get the depth of the quadtree
  unsigned char getDepth() const;

  //! Mask the quadtree
  MaskedQuadtree<T> masked(PolygonMask* m)
  { return MaskedQuadtree<T>(*this, m); }

  friend std::ostream& operator<<<> (std::ostream&, const SmartQuadtree<T>&);

  friend class MaskedQuadtree<T>;
  friend struct const_iterator;
  friend struct iterator;

  friend class Test_SmartQuadtree;

};

template<class T>
struct SmartQuadtree<T>::const_iterator
: std::iterator < std::input_iterator_tag, const T >
{

  const_iterator(
      const typename std::list<SmartQuadtree<T>*>::const_iterator& begin,
      const typename std::list<SmartQuadtree<T>*>::const_iterator& end,
      PolygonMask* mask = NULL);

  const_iterator(const typename SmartQuadtree<T>::iterator& it);

  const_iterator operator++();
  typename SmartQuadtree<T>::const_iterator::reference operator*();
  typename SmartQuadtree<T>::const_iterator::pointer operator->();
  bool operator==(const const_iterator&) const;
  bool operator!=(const const_iterator&) const;

  typename std::vector<const T*>::const_iterator forward_begin();
  typename std::vector<const T*>::const_iterator forward_end();

private:

  typename std::list<SmartQuadtree<T>*>::const_iterator leafIterator, leafEnd;
  typename std::list<T>::const_iterator it, itEnd;
  std::vector<const T*> forward_cells_neighbours;
  typename std::vector<const T*>::const_iterator forward_cells_begin;

  // Current leaf: number of covered summits
  unsigned char aux;
  // NULL if no mask
  PolygonMask* polygonmask;
  // forward_cells_neighbours computed for current cell
  bool neighbours_computed;

  void advanceToNextLeaf();

};


template<class T>
struct SmartQuadtree<T>::iterator
: std::iterator < std::input_iterator_tag, T >
{

  iterator(
      const typename std::list<SmartQuadtree<T>*>::iterator& begin,
      const typename std::list<SmartQuadtree<T>*>::iterator& end,
      PolygonMask* mask = NULL);

  iterator operator++();
  typename SmartQuadtree<T>::iterator::reference operator*();
  typename SmartQuadtree<T>::iterator::pointer operator->();
  bool operator==(const iterator&) const;
  bool operator!=(const iterator&) const;

private:

  typename std::list<SmartQuadtree<T>*>::iterator leafIterator, leafEnd;
  typename std::list<T>::iterator it, itEnd;

  // Elements already parsed
  std::list<T*> already;

  // Current leaf: number of covered summits
  unsigned char aux;
  // NULL if no mask
  PolygonMask* polygonmask;

  void advanceToNextLeaf();

  friend struct SmartQuadtree<T>::const_iterator;
};

template<typename T>
class MaskedQuadtree
{
public:

  MaskedQuadtree(SmartQuadtree<T>& q, PolygonMask* m):
    quadtree(q), polygonmask(m) { }

  typename SmartQuadtree<T>::iterator begin();
  typename SmartQuadtree<T>::iterator end();
  typename SmartQuadtree<T>::const_iterator begin() const;
  typename SmartQuadtree<T>::const_iterator end() const;

private:
  SmartQuadtree<T>& quadtree;
  PolygonMask* polygonmask;

};

#include "quadtree.hpp"

typedef SmartQuadtree<void*> OpaqueSmartQuadtree;
typedef MaskedQuadtree<void*> OpaqueMaskedQuadtree;

#endif // QUADTREE_H
