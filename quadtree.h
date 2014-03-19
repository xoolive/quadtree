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
//   Boundary() {}

  Boundary(float cx, float cy, float dx, float dy) :
    center_x(cx), center_y(cy), dim_x(dx), dim_y(dy),
    x(NULL), y(NULL), limitfct(NULL) {};

  //! Is this point included in the box?
  bool contains(float x, float y);

  //! Is this data included in the box?
  bool contains(void* pt) { return contains(x(pt), y(pt)); }

  //! Max distance between two data in the box
  float norm_l1() { return (dim_x + dim_y); }

  float norm_l2() { return (dim_x * dim_x + dim_y * dim_y); }

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

public:
  // Level differences with the neighbours
  // 0: adjacent quadrant is of same level
  // 1: adjacent quadrant is of larger level (i.e. smaller in size)
  // 2: adjacent quadrant is out of area
  // -1, ...: adjacent quadrant is of smaller level by -n (i.e. larger in size)
  int dn, dw, ds, de;

  // Children nodes
  ExtendedQuadtree *children[4];

private:
  // Data attached to the quadrant
  std::list<void*> points;

  const int capacity;

  // Ancestor
  ExtendedQuadtree* ancestor;

  void updateDN() {
    if (children[2]->samelevel(NORTH)->children[0] != NULL) children[2]->dn = 1;
    if (children[3]->samelevel(NORTH)->children[0] != NULL) children[3]->dn = 1;
  }

  void updateDS() {
    if (children[0]->samelevel(SOUTH)->children[0] != NULL) children[0]->ds = 1;
    if (children[1]->samelevel(SOUTH)->children[0] != NULL) children[1]->ds = 1;
  }

  void updateDE() {
    if (children[1]->samelevel(EAST)->children[0] != NULL) children[1]->de = 1;
    if (children[3]->samelevel(EAST)->children[0] != NULL) children[3]->de = 1;
  }

  void updateDW() {
    if (children[0]->samelevel(WEST)->children[0] != NULL) children[0]->dw = 1;
    if (children[2]->samelevel(WEST)->children[0] != NULL) children[2]->dw = 1;
  }

  bool incrementDS() {
    if (ds < 1) ds += 1;
    if (children[0] == NULL) return false;
    children[0]->incrementDS();
    children[1]->incrementDS();
    return true;
  }
  bool incrementDN() {
    if (dn < 1) dn += 1;
    if (children[0] == NULL) return false;
    children[2]->incrementDN();
    children[3]->incrementDN();
    return true;
  }
  bool incrementDW() {
    if (dw < 1) dw += 1;
    if (children[0] == NULL) return false;
    children[0]->incrementDW();
    children[2]->incrementDW();
    return true;
  }
  bool incrementDE() {
    if (de < 1) de += 1;
    if (children[0] == NULL) return false;
    children[1]->incrementDE();
    children[3]->incrementDE();
    return true;
  }

public:

  //! Constructor
  ExtendedQuadtree(float center_x, float center_y, float dim_x, float dim_y,
                   int capacity) :
    b(center_x, center_y, dim_x, dim_y), location(0), level(0),
    dn(2), dw(2), ds(2), de(2), capacity(capacity), ancestor(this)
  {
    for (int i = 0; i<4; ++i) children[i] = NULL;
  }

  //! Constructor of a child quadtree
  // SW -> 0, SE -> 1, NW -> 2, NE -> 3
  ExtendedQuadtree(const ExtendedQuadtree&, int);

  //! Find same level neighbour in determined direction
  ExtendedQuadtree* samelevel(enum Direction) const;

  //! Insert one piece of data to the quadrant
  bool insert(void*);

  //! Returns the subquadrant pointed by location code
  ExtendedQuadtree* getQuadrant(std::size_t location, std::size_t level) const;

  //! Returns the data embedded to current quadrant
  const std::list<void*>& getPoints() const { return points; }

  //! Sets the getters in the Boundary
  void setXYFcts(float (*x)(void*), float (*y)(void*)) { b.x = x; b.y = y; }

  void setLimitation(bool (*limitation)(Boundary*)) { b.limitfct = limitation; }

  //! Get the location
  std::size_t getLocation() const { return location; }

  //! Get the level
  std::size_t getLevel() const { return level; }

  friend std::ostream& operator<<(std::ostream&, const ExtendedQuadtree&);

};

#endif // QUADTREE_H
