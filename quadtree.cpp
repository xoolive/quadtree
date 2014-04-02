#include "quadtree.h"

#include <cfloat> // FLT_EPSILON
#include <stack>
#include <vector>

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
  /*static*/ std::stack<unsigned char> stack;

  for (size_t i = 0; i < l; i++, location >>= 2)
    stack.push(location & 3);
  for (size_t i = 0; i < l; i++)
  {
    desc = quadrant->children[stack.top()];
    if (NULL == desc) {
      // never much, but avoid memory leaks (see static)
      // while (!stack.empty()) stack.pop();
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

  assert(b.limitfct != NULL);

  for (int i = 0; i<4; ++i) children[i] = NULL;

  location = (e.location << 2) + direction;
  level    = e.level + 1;
  ancestor = e.ancestor;

  // NSWE updates

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

  // Diagonal updates

  if (direction == 0)
  {
    delta[NORTHEAST] = 0;
    delta[NORTHWEST] = 3;
    delta[SOUTHWEST] = (e.delta[SOUTHWEST]>1 ? e.delta[SOUTHWEST] :
                        e.delta[SOUTHWEST]-1);;
    delta[SOUTHEAST] = 3;
  }
  else if (direction == 1)
  {
    delta[NORTHWEST] = 0;
    delta[SOUTHWEST] = 3;
    delta[SOUTHEAST] = (e.delta[SOUTHEAST]>1 ? e.delta[SOUTHEAST] :
                        e.delta[SOUTHEAST]-1);
    delta[NORTHEAST] = 3;
  }
  else if (direction == 2)
  {
    delta[SOUTHEAST] = 0;
    delta[NORTHEAST] = 3;
    delta[NORTHWEST] = (e.delta[NORTHWEST]>1 ? e.delta[NORTHWEST] :
                        e.delta[NORTHWEST]-1);
    delta[SOUTHWEST] = 3;
  }
  else // direction == 3
  {
    delta[SOUTHWEST] = 0;
    delta[SOUTHEAST] = 3;
    delta[NORTHEAST] = (e.delta[NORTHEAST]>1 ? e.delta[NORTHEAST] :
                        e.delta[NORTHEAST]-1);
    delta[NORTHWEST] = 3;
  }

  b.dim_x  = e.b.dim_x / 2.;
  b.dim_y  = e.b.dim_y / 2.;

  b.limit = b.limitation();

}

bool ExtendedQuadtree::insert(void* pt)
{
  if (!b.contains(pt)) return false;

  // It is OK to go over capacity if a test "limitation" on b is verified
  if (b.limit || ((NULL == children[0]) && (points.size() < capacity)))
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

void ExtendedQuadtree::updateDiagonal(unsigned int diagdir,
                                      unsigned int dir, int d)
{
  if (children[0] == NULL)
  {
    assert(delta[diagdir] == 3);
    delta[diagdir] = d;
    this->samelevel(diagdir)->delta[ (diagdir+4)&7 ] = (0==d?d:1);
    return ;
  }

  // TODO make it easier to read?

  if (dir == WEST)
  {
    if (diagdir == NORTHWEST)
      children[2]->updateDiagonal(diagdir, dir, d-1);
    else
      children[0]->updateDiagonal(diagdir, dir, d-1);
  }
  if (dir == SOUTH)
  {
    if (diagdir == SOUTHEAST)
      children[1]->updateDiagonal(diagdir, dir, d-1);
    else
      children[0]->updateDiagonal(diagdir, dir, d-1);
  }
  if (dir == EAST)
  {
    if (diagdir == NORTHEAST)
      children[3]->updateDiagonal(diagdir, dir, d-1);
    else
      children[1]->updateDiagonal(diagdir, dir, d-1);
  }
  if (dir == NORTH)
  {
    if (diagdir == NORTHEAST)
      children[3]->updateDiagonal(diagdir, dir, d-1);
    else
      children[2]->updateDiagonal(diagdir, dir, d-1);
  }

}

//! Increments the delta in direction dir
//! Returns true if you have children
bool ExtendedQuadtree::incrementDelta(unsigned int dir, bool flag)
{
  if (children[0] == NULL)
  {
    if (delta[dir] < 1) delta[dir] += 1;
    return false;
  }

  if (flag)
  {
    if (dir == WEST)
    {
      children[0]->updateDiagonal(NORTHWEST, dir, 0);
      children[2]->updateDiagonal(SOUTHWEST, dir, 0);
    }
    if (dir == SOUTH)
    {
      children[0]->updateDiagonal(SOUTHEAST, dir, 0);
      children[1]->updateDiagonal(SOUTHWEST, dir, 0);
    }
    if (dir == EAST)
    {
      children[1]->updateDiagonal(NORTHEAST, dir, 0);
      children[3]->updateDiagonal(SOUTHEAST, dir, 0);
    }
    if (dir == NORTH)
    {
      children[2]->updateDiagonal(NORTHEAST, dir, 0);
      children[3]->updateDiagonal(NORTHWEST, dir, 0);
    }
  }


  if ( dir < 3 ) // NORTHEAST
    children[3]->incrementDelta(dir, false);
  if ( ((dir + 6) & 7) < 3 ) // NORTHWEST
    children[2]->incrementDelta(dir, false);
  if ( ((dir + 4) & 7) < 3 ) // SOUTHWEST
    children[0]->incrementDelta(dir, false);
  if ( ((dir + 2) & 7) < 3 ) // SOUTHEAST
    children[1]->incrementDelta(dir, false);
  return true;
}

//! Updates the delta in direction dir if neighbour of same level has
//! children nodes
void ExtendedQuadtree::updateDelta(unsigned int dir)
{
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

void ExtendedQuadtree::iterate(bool (*apply)(void*))
{
  if (children[0] != NULL)
  {
    children[0]->iterate(apply);
    children[1]->iterate(apply);
    children[2]->iterate(apply);
    children[3]->iterate(apply);
  }

  std::list<void*>::iterator it = points.begin(), ie = points.end();
  for ( ; it != ie; ++it)
    if (apply(*it))
      if (!b.contains(*it))
      {
        // Computing the proper neighbour is as fast as finding it from the
        // ancestor node...
        ancestor->insert(*it);
        it = points.erase(it);
      }
}

void ExtendedQuadtree::iterate(void (*apply)(void*, void*))
{
  if (children[0] != NULL)
  {
    children[0]->iterate(apply);
    children[1]->iterate(apply);
    children[2]->iterate(apply);
    children[3]->iterate(apply);
  }

  std::vector<void*> neighbours;
  std::list<void*>::const_iterator it = points.begin(), ie = points.end();
  for ( ; it != ie; ++it) { neighbours.push_back((*it)); }

  ExtendedQuadtree* nb;
  for (size_t i = 0; i < 4; ++i)
    if ( (i < 4 && delta[i] < 1) || delta[i] < 0)
    {
      nb = samelevel(i);
      it = nb->getPoints().begin(), ie = nb->getPoints().end();
      for ( ; it != ie; ++it) { neighbours.push_back((*it)); }
    }

  for (int i = 0; i < points.size(); ++i)
    for (int j = i+1; j < neighbours.size(); ++j)
      apply(neighbours[i], neighbours[j]);

}

std::size_t ExtendedQuadtree::getDataSize() const
{
  std::size_t size = 0, tmp;
  if (children[0] != NULL)
  {
    tmp = children[0]->getDataSize();
    if (tmp > size) size = tmp;
    tmp = children[1]->getDataSize();
    if (tmp > size) size = tmp;
    tmp = children[2]->getDataSize();
    if (tmp > size) size = tmp;
    tmp = children[3]->getDataSize();
    if (tmp > size) size = tmp;
    return size;
  }

  return points.size();
}

std::size_t ExtendedQuadtree::getDepth() const
{
  std::size_t depth = 0, tmp;
  if (children[0] != NULL)
  {
    tmp = children[0]->getDepth();
    if (tmp > depth) depth = tmp;
    tmp = children[1]->getDepth();
    if (tmp > depth) depth = tmp;
    tmp = children[2]->getDepth();
    if (tmp > depth) depth = tmp;
    tmp = children[3]->getDepth();
    if (tmp > depth) depth = tmp;
    return depth + 1;
  }

  return 0;
}

