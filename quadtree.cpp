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
ExtendedQuadtree::samelevel(unsigned int dir) const
{
  if (delta[dir] == 2) return NULL;
  unsigned int newloc = Neighbour::samelevel(location, dir, level);
  return getQuadrant(newloc, level);
}

ExtendedQuadtree*
ExtendedQuadtree::getQuadrant(std::size_t location, std::size_t l) const
{
  ExtendedQuadtree *quadrant = ancestor, *desc;
  static std::stack<unsigned char> stack;

  for (size_t i = 0; i < l; i++, location >>= 2)
    stack.push(location & 3);
  for (size_t i = 0; i < l; i++)
  {
    desc = quadrant->children[stack.top()];
    if (NULL == desc) {
      // never much, but avoid memory leaks (see static)
      while (!stack.empty()) stack.pop();
      return quadrant;
    }
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
    delta[SOUTH] = 0;
    delta[NORTH] = (e.delta[NORTH]==2 ? 2 : e.delta[NORTH]-1);
  }
  else // south
  {
    b.center_y = e.b.center_y - e.b.dim_y / 2.;
    delta[NORTH] = 0;
    delta[SOUTH] = (e.delta[SOUTH]==2 ? 2 : e.delta[SOUTH]-1);
  }
  if (((direction) & 1) == 0) // west
  {
    b.center_x = e.b.center_x - e.b.dim_x / 2.;
    delta[EAST] = 0;
    delta[WEST] = (e.delta[WEST]==2 ? 2 : e.delta[WEST]-1);
  }
  else // east
  {
    b.center_x = e.b.center_x + e.b.dim_x / 2.;
    delta[WEST] = 0;
    delta[EAST] = (e.delta[EAST]==2 ? 2 : e.delta[EAST]-1);
  }

  if (direction == 0)
  {
    delta[NORTHEAST] = 0;
    delta[NORTHWEST] = delta[WEST];
    delta[SOUTHWEST] = (e.delta[SOUTHWEST]==2 ? 2 : e.delta[SOUTHWEST]-1);;
    delta[SOUTHEAST] = delta[SOUTH];
  }
  else if (direction == 1)
  {
    delta[NORTHWEST] = 0;
    delta[SOUTHWEST] = delta[SOUTH];
    delta[SOUTHEAST] = (e.delta[SOUTHEAST]==2 ? 2 : e.delta[SOUTHEAST]-1);
    delta[NORTHEAST] = delta[EAST];
  }
  else if (direction == 2)
  {
    delta[SOUTHEAST] = 0;
    delta[NORTHEAST] = delta[NORTH];
    delta[NORTHWEST] = (e.delta[NORTHWEST]==2 ? 2 : e.delta[NORTHWEST]-1);
    delta[SOUTHWEST] = delta[WEST];
  }
  else // direction == 3
  {
    delta[SOUTHWEST] = 0;
    delta[SOUTHEAST] = delta[EAST];
    delta[NORTHEAST] = (e.delta[NORTHEAST]==2 ? 2 : e.delta[NORTHEAST]-1);
    delta[NORTHWEST] = delta[NORTH];
  }

  b.dim_x  = e.b.dim_x / 2.;
  b.dim_y  = e.b.dim_y / 2.;

}

bool ExtendedQuadtree::insert(void* pt)
{
  if (!b.contains(pt)) return false;

  // It is OK to go over capacity if a test "limitation" on b is verified
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

    std::cout << "increment all for 0x" << std::hex << location << " " <<
      std::dec << level <<  std::endl;
    // Update neighbour info
    // TODO find a solution for not parsing all directions!
    for (unsigned int i = 0; i < 8; ++i)
      if (this->delta[i] < 2)
        if (this->samelevel(i)->incrementDelta((i+4) & 7))
          updateDelta(i);

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

//! Increments the delta in direction dir
//! Returns true if you have children
bool ExtendedQuadtree::incrementDelta(unsigned int dir)
{
  std::cout << "increment delta for 0x" << std::hex << location << " " <<
    std::dec << level << " " << dir << std::endl;
  if (children[0] == NULL)
  {
    if (delta[dir] < 1) delta[dir] += 1;
    return false;
  }
  if ( dir < 3 ) // NORTHEAST
    children[3]->incrementDelta(dir);
  if ( ((dir + 6) & 7) < 3 ) // NORTHWEST
    children[2]->incrementDelta(dir);
  if ( ((dir + 4) & 7) < 3 ) // SOUTHWEST
    children[0]->incrementDelta(dir);
  if ( ((dir + 2) & 7) < 3 ) // SOUTHEAST
    children[1]->incrementDelta(dir);
  return true;
}

//! Updates the delta in direction dir if neighbour of same level has
//! children nodes
void ExtendedQuadtree::updateDelta(unsigned int dir)
{
  std::cout << "update delta for 0x" << std::hex <<location << " " <<std::dec << level << " " << dir <<std::endl;
  if ( dir < 3 ) // NORTHEAST corner
    if (children[3]->samelevel(dir)->children[0] != NULL)
      children[3]->delta[dir] = 1;
  if ( ((dir + 6) & 7) < 3 ) // NORTHWEST corner
    if (children[2]->samelevel(dir)->children[0] != NULL)
      children[2]->delta[dir] = 1;
  if ( ((dir + 4) & 7) < 3 ) // SOUTHWEST corner
    if (children[0]->samelevel(dir)->children[0] != NULL)
      children[0]->delta[dir] = 1;
  if ( ((dir + 2) & 7) < 3 ) // SOUTHEAST corner
    if (children[1]->samelevel(dir)->children[0] != NULL)
      children[1]->delta[dir] = 1;
}

