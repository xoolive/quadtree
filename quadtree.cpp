#include "quadtree.h"

#include <cfloat> // FLT_EPSILON
#include <stack>
#include <vector>

const unsigned char ExtendedQuadtree::diags[] =
{ SOUTHWEST, SOUTHEAST, NORTHWEST, NORTHEAST };

bool Boundary::contains(float x, float y)
{
  return ((x < center_x + dim_x + FLT_EPSILON) &&
          (x > center_x - dim_x - FLT_EPSILON) &&
          (y < center_y + dim_y + FLT_EPSILON) &&
          (y > center_y - dim_y - FLT_EPSILON));
}

ExtendedQuadtree*
ExtendedQuadtree::samelevel(unsigned char dir) const
{
  if (delta[dir] == 2) return NULL;
  unsigned int newloc = Neighbour::samelevel(location, dir, level);
  return getQuadrant(newloc, level);
}

ExtendedQuadtree*
ExtendedQuadtree::getQuadrant(unsigned long location,
                              unsigned short depth) const
{
  assert(depth < 2048);
  ExtendedQuadtree *quadrant = ancestor, *desc;
  static unsigned char stack[2048];
  short istack = depth - 1;

  for (unsigned short i = 0; i < depth; i++, location >>= 2)
    stack[i] = location & 3;
  for (unsigned short i = 0; i < depth; i++)
  {
    desc = quadrant->children[stack[istack]];
    if (NULL == desc)
      return quadrant;
    quadrant = desc;
    istack --;
  }

  return quadrant;
}

ExtendedQuadtree::ExtendedQuadtree(const ExtendedQuadtree& e,
                                   unsigned char subdivision)
  : b(e.b), capacity(e.capacity)
{

  assert(b.limitfct != NULL);

  for (int i = 0; i<4; ++i) children[i] = NULL;

  location = (e.location << 2) + subdivision;
  level    = e.level + 1;
  ancestor = e.ancestor;

  if (subdivision > 1) // north
    b.center_y = e.b.center_y + e.b.dim_y / 2.;
  else // south
    b.center_y = e.b.center_y - e.b.dim_y / 2.;
  if (((subdivision) & 1) == 0) // west
    b.center_x = e.b.center_x - e.b.dim_x / 2.;
  else // east
    b.center_x = e.b.center_x + e.b.dim_x / 2.;

  // Updating delta
  const unsigned char diag = diags[subdivision];
  delta[diag] = (e.delta[diag]> 1 ? e.delta[diag] : e.delta[diag]-1);
  delta[(diag+1)&7] = (e.delta[(diag+1)&7]==2 ? 2 : e.delta[(diag+1)&7]-1);
  delta[(diag+2)&7] = 3;
  delta[(diag+3)&7] = 0;
  delta[(diag+4)&7] = 0;
  delta[(diag+5)&7] = 0;
  delta[(diag+6)&7] = 3;
  delta[(diag+7)&7] = (e.delta[(diag+7)&7]==2 ? 2 : e.delta[(diag+7)&7]-1);

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

void ExtendedQuadtree::updateDiagonal(unsigned char diagdir,
                                      unsigned char dir, int d)
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
bool ExtendedQuadtree::incrementDelta(unsigned char dir, bool flag)
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
void ExtendedQuadtree::updateDelta(unsigned char dir)
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
    return;
  }

  std::list<void*>::iterator it = points.begin(), ie = points.end();
  while ( it != ie )
    if (apply(*it))
      if (!b.contains(*it))
      {
        // Computing the proper neighbour is as fast as finding it from the
        // ancestor node...
        ancestor->insert(*it);
        it = points.erase(it);
      }
      else ++it;
    else ++it;
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

  ExtendedQuadtree* nb;
  for (size_t i = 0; i < 4; ++i)
    if (delta[i] < 1) {
      nb = samelevel(i);
      neighbours.insert(neighbours.end(),
                        nb->getPoints().begin(), nb->getPoints().end());
    }
  for (size_t i = 4; i < 8; ++i)
    if (delta[i] < 0)
    {
      nb = samelevel(i);
      neighbours.insert(neighbours.end(),
                        nb->getPoints().begin(), nb->getPoints().end());
    }

  std::list<void*>::const_iterator it = points.begin(), ie = points.end();
  for ( ; it != ie ; ++it)
  {
    std::list<void*>::const_iterator jt = it; ++jt;
    for ( ; jt != ie; ++jt)
      apply(*it, *jt);
    for (int j = 0; j < neighbours.size(); ++j)
      apply(*it, neighbours[j]);
  }

}

void ExtendedQuadtree::iterateby4(void (*apply)(void*, void*),
                                  void (*applyby4)(void*, void**))
{
  if (children[0] != NULL)
  {
    children[0]->iterateby4(apply, applyby4);
    children[1]->iterateby4(apply, applyby4);
    children[2]->iterateby4(apply, applyby4);
    children[3]->iterateby4(apply, applyby4);
  }

  std::vector<void*> neighbours;
  std::list<void*>::const_iterator it = points.begin(), ie = points.end();

  ExtendedQuadtree* nb;
  for (size_t i = 0; i < 4; ++i)
    if (delta[i] < 1) {
      nb = samelevel(i);
      it = nb->getPoints().begin(), ie = nb->getPoints().end();
      neighbours.insert(neighbours.end(),it,ie);
    }
  for (size_t i = 4; i < 8; ++i)
    if (delta[i] < 0)
    {
      nb = samelevel(i);
      it = nb->getPoints().begin(), ie = nb->getPoints().end();
      neighbours.insert(neighbours.end(),it,ie);
    }

    for (std::list<void*>::iterator itL = points.begin(); itL != points.end();
         itL++) {
      std::list<void*>::iterator jtL = itL;
      jtL++;
      for ( ; jtL != points.end(); ++jtL)
        apply(*itL, *jtL);
      int j=0;
      for ( ; j < int(neighbours.size())-3; j+=4)
        applyby4(*itL, &(neighbours[j]));
      for( ; j < neighbours.size(); j++)
        apply(*itL, neighbours[j]);
    }
}


unsigned long ExtendedQuadtree::getDataSize() const
{
  unsigned long size = 0, tmp;
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

unsigned char ExtendedQuadtree::getDepth() const
{
  unsigned char depth = 0, tmp;
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

