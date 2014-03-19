#include "quadtree.h"

#include <cfloat> // FLT_EPSILON
#include <stack>

bool Boundary::contains(float x, float y)
{
  return ((x < center_x + dim_x + FLT_EPSILON) &&
          (x > center_x - dim_x - FLT_EPSILON) &&
          (y < center_y + dim_y + FLT_EPSILON) &&
          (y > center_y - dim_y - FLT_EPSILON));
}

ExtendedQuadtree*
ExtendedQuadtree::samelevel(enum Direction dir) const
{
  switch (dir) {
  case NORTH:
  case NORTHEAST:
  case NORTHWEST:
    if (dn == 2) return NULL;
    break;
  case EAST:
  case SOUTHEAST:
    if (de == 2) return NULL;
    break;
  case SOUTH:
  case SOUTHWEST:
    if (ds == 2) return NULL;
    break;
  default:
    if (dw == 2) return NULL;
  }

  unsigned int newloc = Neighbour::samelevel(location, dir, level);

  return getQuadrant(newloc, level);
}

ExtendedQuadtree*
ExtendedQuadtree::getQuadrant(std::size_t location, std::size_t l) const
{
  ExtendedQuadtree *quadrant = ancestor, *desc;
  static std::stack<unsigned char> stack;

  for (size_t i = 0; i < l; i++, location >>=2)
    stack.push(location & 3);
  for (size_t i = 0; i < l; i++)
  {
    desc = quadrant->children[stack.top()];
    if (NULL == desc) return quadrant;
    quadrant = desc;
    stack.pop();
  }

  return quadrant;
}

ExtendedQuadtree::ExtendedQuadtree(const ExtendedQuadtree& e, int direction)
  : b(e.b), capacity(e.capacity)
{
  for (int i = 0; i<4; ++i) children[i] = NULL;
  location = (e.location << 2) + direction;
  level    = e.level + 1;
  ancestor = e.ancestor;

  if (direction > 1) // north
  {
    b.center_y = e.b.center_y + e.b.dim_y / 2.;
    ds = 0;
    dn = (e.dn==2?2:e.dn-1);
  }
  else
  {
    b.center_y = e.b.center_y - e.b.dim_y / 2.;
    dn = 0;
    ds = (e.ds==2?2:e.ds-1);
  }
  if (((direction) & 1) == 0) // west
  {
    b.center_x = e.b.center_x - e.b.dim_x / 2.;
    de = 0;
    dw = (e.dw==2?2:e.dw-1);
  }
  else
  {
    b.center_x = e.b.center_x + e.b.dim_x / 2.;
    dw = 0;
    de = (e.de==2?2:e.de-1);
  }

  b.dim_x  = e.b.dim_x / 2.;
  b.dim_y  = e.b.dim_y / 2.;

}

bool ExtendedQuadtree::insert(void* pt)
{
  if (!b.contains(pt)) return false;

  // It is OK to go over capacity if a test on b is verified
  if (b.limitation() || ((NULL == children[0]) && (points.size() < capacity)))
  {
    points.push_back(pt);
    return true;
  }

  if (NULL == children[0])
  {
    children[0] = new ExtendedQuadtree(*this, 0);
    children[1] = new ExtendedQuadtree(*this, 1);
    children[2] = new ExtendedQuadtree(*this, 2);
    children[3] = new ExtendedQuadtree(*this, 3);

    // Update neighbour info
    if (this->dn < 2) 
      if (this->samelevel(NORTH)->incrementDS()) updateDN();
    if (this->ds < 2)
      if (this->samelevel(SOUTH)->incrementDN()) updateDS();
    if (this->de < 2)
      if (this->samelevel(EAST)->incrementDW()) updateDE();
    if (this->dw < 2)
      if (this->samelevel(WEST)->incrementDE()) updateDW();

    // Forward data to children
    std::list<void*>::iterator it = points.begin(), ie = points.end();
    for ( ; it != ie; ++it) this->insert(*it);
    points.clear();
  }

  if (children[0]->insert(pt)) return true;
  if (children[1]->insert(pt)) return true;
  if (children[2]->insert(pt)) return true;
  if (children[3]->insert(pt)) return true;

  return false;

}
