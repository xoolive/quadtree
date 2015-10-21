/*
 * Interface for a smart version of quadtrees specialised for tracking moving
 * objects. Iterators and their const counterparts let you process objects and
 * ajust the quadtree structure when points are moved.
 *
 * Xavier Olive, 28 nov. 2014
 */

#include <vector>
#include <algorithm>

using std::list;

template<typename T>
const unsigned char SmartQuadtree<T>::diags[] =
{ SOUTHWEST, SOUTHEAST, NORTHWEST, NORTHEAST };

template<typename T>
SmartQuadtree<T>* SmartQuadtree<T>::samelevel(unsigned char dir) const
{
  if (delta[dir] == 2) return NULL;
  unsigned int newloc = Neighbour::samelevel(location, dir, level);
  return getQuadrant(newloc, level);
}

template<typename T>
SmartQuadtree<T>* SmartQuadtree<T>::getQuadrant(unsigned long location,
                                                unsigned short depth) const
{
  assert(depth < 2048);
  SmartQuadtree *quadrant = ancestor, *desc;
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

template<typename T>
SmartQuadtree<T>::SmartQuadtree(const SmartQuadtree<T>& e,
                                unsigned char subdivision,
                                typename list<SmartQuadtree<T>*>::iterator& w)
: b(e.b), capacity(e.capacity)
{

  for (int i = 0; i<4; ++i) children[i] = NULL;

  location = (e.location << 2) + subdivision;
  level    = e.level + 1;
  ancestor = e.ancestor;

  w = ancestor->leaves.insert(w, this);

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
  delta[(diag+1) & 7] = (e.delta[(diag+1)&7]==2 ? 2 : e.delta[(diag+1)&7]-1);
  delta[(diag+2) & 7] = 3;
  delta[(diag+3) & 7] = 0;
  delta[(diag+4) & 7] = 0;
  delta[(diag+5) & 7] = 0;
  delta[(diag+6) & 7] = 3;
  delta[(diag+7) & 7] = (e.delta[(diag+7)&7]==2 ? 2 : e.delta[(diag+7)&7]-1);

  b.dim_x  = e.b.dim_x / 2.;
  b.dim_y  = e.b.dim_y / 2.;

  b.limit = BoundaryLimit<Boundary>::limitation(b);

}

template<typename T>
SmartQuadtree<T>::~SmartQuadtree()
{
  if (NULL == children[0]) return;
  delete children[0]; delete children[1];
  delete children[2]; delete children[3];
}

template<typename T>
typename SmartQuadtree<T>::const_iterator SmartQuadtree<T>::begin() const
{ return SmartQuadtree<T>::const_iterator(leaves.begin(), leaves.end()); }

template<typename T>
typename SmartQuadtree<T>::const_iterator SmartQuadtree<T>::end() const
{ return SmartQuadtree<T>::const_iterator(leaves.end(), leaves.end()); }

template<typename T>
typename SmartQuadtree<T>::iterator SmartQuadtree<T>::begin()
{ return SmartQuadtree<T>::iterator(leaves.begin(), leaves.end()); }

template<typename T>
typename SmartQuadtree<T>::iterator SmartQuadtree<T>::end()
{ return SmartQuadtree<T>::iterator(leaves.end(), leaves.end()); }

template<typename T>
bool SmartQuadtree<T>::insert(T pt)
{
  if (!b.contains(pt)) return false;

  // It is OK to go over capacity if a test "limitation" on b is verified
  if (b.limit || ((NULL == children[0]) && (points.size() < capacity)))
  {
    points.push_back(pt);
    ancestor->where[&points.back()] = this;
    return true;
  }

  if (NULL == children[0])
  {
    typename list<SmartQuadtree<T>*>::iterator w =
      std::find(ancestor->leaves.begin(), ancestor->leaves.end(), this);

    assert(w != ancestor->leaves.end());
    w = ancestor->leaves.erase(w);

    children[0] = new SmartQuadtree(*this, 0, w);
    children[1] = new SmartQuadtree(*this, 1, w);
    children[2] = new SmartQuadtree(*this, 2, w);
    children[3] = new SmartQuadtree(*this, 3, w);

    // Update neighbour info
    for (unsigned int i = 0; i < 8; ++i)
      if (this->delta[i] < 2)
        if (this->samelevel(i)->incrementDelta((i+4) & 7))
          updateDelta(i);

    // Forward data to children
    //std::for_each (points.begin(), points.end(),
    //               [=](void* p) { this->insert(p); } );
    for (typename list<T>::iterator it = points.begin();
        it != points.end(); ++it)
      this->insert(*it);
    points.clear();
  }

  if (children[0]->insert(pt)) return true;
  if (children[1]->insert(pt)) return true;
  if (children[2]->insert(pt)) return true;
  if (children[3]->insert(pt)) return true;

  return false;

}

template<typename T>
void SmartQuadtree<T>::updateDiagonal(unsigned char diagdir,
                                      unsigned char dir, int d)
{
  if (children[0] == NULL)
  {
    assert(delta[diagdir] == 3);
    delta[diagdir] = d;
    this->samelevel(diagdir)->delta[ (diagdir+4)&7 ] = (0==d?d:1);
    return ;
  }

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
template<typename T>
bool SmartQuadtree<T>::incrementDelta(unsigned char dir, bool flag)
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
template<typename T>
void SmartQuadtree<T>::updateDelta(unsigned char dir)
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
/*
template<typename T>
void SmartQuadtree<T>::iterate(const PolygonMask& m, bool (*apply)(T&))
{
  PolygonMask clip = m.clip(b);
  if (clip.getSize() < 3) return;

  int nb = b.coveredByPolygon(clip);

  // the quadtree is inside the polygon
  if (nb == 4) {
      //return iterate(apply);
      SmartQuadtree<T>::iterator it = begin();
      for (; it != end(); ++it)
          apply(*it);
  }

  if (children[0] != NULL)
  {
    children[0]->iterate(clip, apply);
    children[1]->iterate(clip, apply);
    children[2]->iterate(clip, apply);
    children[3]->iterate(clip, apply);
    return;
  }

  typename list<T>::iterator it = points.begin(), ie = points.end();
  while ( it != ie )
    if (m.pointInPolygon(b.x(&(*it)), b.y(&(*it))))
      if (apply(*it))
        if (!b.contains(&(*it)))
        {
          // Computing the proper neighbour is as fast as finding it from the
          // ancestor node...
          ancestor->insert(*it);
          it = points.erase(it);
        }
        else ++it;
      else ++it;
    else ++it;

}
*/
/*template<typename T>
bool SmartQuadtree<T>::updateData(T& p)
{
  SmartQuadtree* e = ancestor->where[&p];
  assert (e != NULL);
  assert (e->points.end() != std::find(e->points.begin(), e->points.end(), p));

  if (e->contains(p)) return false;
  e->points.remove(p);
  ancestor->insert(p);
  return true;
}
*/
template<typename T>
void SmartQuadtree<T>::removeData(T& p)
{
  SmartQuadtree* e = ancestor->where[&p];
  assert (e != NULL);
  assert (e->points.end() != std::find(e->points.begin(), e->points.end(), p));

  e->points.remove(p);
  ancestor->where.erase(p);
}

template<typename T>
SmartQuadtree<T>::const_iterator::const_iterator(
    const typename list<SmartQuadtree<T>*>::const_iterator& begin,
    const typename list<SmartQuadtree<T>*>::const_iterator& end,
    PolygonMask* mask) : polygonmask(mask)
{
  leafIterator = begin;
  leafEnd = end;
  aux = 4; // Default case: polygonmask is not set
  if (begin != end)
  {
    it = (*leafIterator)->points.begin();
    itEnd = (*leafIterator)->points.end();
    if (polygonmask != NULL)
    {
      // In case the first leaf is out of the polygon
      PolygonMask clip = polygonmask->clip((*leafIterator)->b);
      aux = (*leafIterator)->b.coveredByPolygon(clip);
    }
    // In case it = itEnd
    advanceToNextLeaf();
    assert(leafIterator != leafEnd ? it != itEnd : true);
    if (aux < 4)
      // This only happens if polygonmask is set
      while (!polygonmask->pointInPolygon((*leafIterator)->b.x(&(*it)),
                                          (*leafIterator)->b.y(&(*it))))
      {
        ++it;
        advanceToNextLeaf();
      }
    assert(leafIterator != leafEnd ? it != itEnd : true);
  }
}

template<typename T>
SmartQuadtree<T>::const_iterator::const_iterator(
    const typename SmartQuadtree<T>::iterator& a)
{
  leafIterator = a.leafIterator;
  leafEnd = a.leafEnd;
  it = a.it;
  itEnd = a.itEnd;
  polygonmask = a.polygonmask;
  aux = a.aux;
}

template<typename T>
void SmartQuadtree<T>::const_iterator::advanceToNextLeaf()
{
  if (it == itEnd)
    do
    {
      ++leafIterator;
      if (leafIterator == leafEnd) return;
      if (polygonmask != NULL)
      {
        PolygonMask clip = polygonmask->clip((*leafIterator)->b);
        if (clip.getSize() < 3) continue;
        aux = (*leafIterator)->b.coveredByPolygon(clip);
      }
      it = (*leafIterator)->points.begin();
      itEnd = (*leafIterator)->points.end();
    } while (it == itEnd);
}

template<typename T>
typename SmartQuadtree<T>::const_iterator
SmartQuadtree<T>::const_iterator::operator++()
{
  if (leafIterator == leafEnd) return *this;
  assert (it != itEnd);
  ++it;
  advanceToNextLeaf();
  if (aux < 4)
  {
    // If a polygonmask is set, we want to ensure than (*it) is inside
    assert(polygonmask != NULL);
    while ((leafIterator != leafEnd) &&
           !polygonmask->pointInPolygon(BoundaryXY<T>::getX(*it),
                                        BoundaryXY<T>::getY(*it))
          )
    {
      ++it;
      advanceToNextLeaf();
    }
  }
  assert(leafIterator != leafEnd ? it != itEnd : true);
  return *this;
}

template<typename T>
typename SmartQuadtree<T>::const_iterator::reference
SmartQuadtree<T>::const_iterator::operator*()
{ return it.operator*(); }

template<typename T>
typename SmartQuadtree<T>::const_iterator::pointer
SmartQuadtree<T>::const_iterator::operator->()
{ return it.operator->(); }

template<typename T> bool
SmartQuadtree<T>::const_iterator::operator==(
    const typename SmartQuadtree<T>::const_iterator& rhs) const
{
  if (it == itEnd) return leafIterator == rhs.leafIterator;
  return (it == rhs.it) && (leafIterator == rhs.leafIterator);
}

template<typename T> bool
SmartQuadtree<T>::const_iterator::operator!=(
    const typename SmartQuadtree<T>::const_iterator& rhs) const
{ return !(*this == rhs); }


template<typename T>
SmartQuadtree<T>::iterator::iterator(
    const typename list<SmartQuadtree<T>*>::iterator& begin,
    const typename list<SmartQuadtree<T>*>::iterator& end,
    PolygonMask* mask) : polygonmask(mask)
{
  leafIterator = begin;
  leafEnd = end;
  aux = 4; // Default case: polygonmask is not set
  if (begin != end)
  {
    it = (*leafIterator)->points.begin();
    itEnd = (*leafIterator)->points.end();
    if (polygonmask != NULL)
    {
      // In case the first leaf is out of the polygon
      PolygonMask clip = polygonmask->clip((*leafIterator)->b);
      aux = (*leafIterator)->b.coveredByPolygon(clip);
    }
    // In case it = itEnd
    advanceToNextLeaf();
    assert(leafIterator != leafEnd ? it != itEnd : true);
    if (aux < 4)
      // This only happens if polygonmask is set
      while (!polygonmask->pointInPolygon(BoundaryXY<T>::getX(*it),
                                          BoundaryXY<T>::getY(*it)))
      {
        ++it;
        advanceToNextLeaf();
      }
    assert(leafIterator != leafEnd ? it != itEnd : true);
  }
}

template<typename T>
void SmartQuadtree<T>::iterator::advanceToNextLeaf()
{
  if (it == itEnd)
    do
    {
      ++leafIterator;
      if (leafIterator == leafEnd) return;
      if (polygonmask != NULL)
      {
        PolygonMask clip = polygonmask->clip((*leafIterator)->b);
        if (clip.getSize() < 3) continue;
        aux = (*leafIterator)->b.coveredByPolygon(clip);
      }
      it = (*leafIterator)->points.begin();
      itEnd = (*leafIterator)->points.end();
    } while (it == itEnd);
}


template<typename T>
typename SmartQuadtree<T>::iterator
SmartQuadtree<T>::iterator::operator++()
{
  if (leafIterator == leafEnd) return *this;
  assert (it != itEnd);
  if (!(*leafIterator)->b.contains(*it))
  {
    // Computing the proper neighbour is probably slower than finding it
    // from the ancestor node...
    SmartQuadtree<T>* previous = (*leafIterator)->ancestor->where[&(*it)];
    assert (previous != NULL);
    (*leafIterator)->ancestor->insert(*it);
    SmartQuadtree<T>* current = (*leafIterator)->ancestor->where[&(*it)];
    assert (current != NULL);
    if (current->location > previous->location)
      already.push_back(&(*it));
    it = (*leafIterator)->points.erase(it);
  }
  else
    ++it;

  advanceToNextLeaf();

  // Don't parse elements that are already parsed
  if (aux == 4)
  while (already.end() != std::find(already.begin(), already.end(), &(*it)))
  {
    ++it;
    advanceToNextLeaf();
  }
  if (aux < 4)
  {
    // If a polygonmask is set, we want to ensure than (*it) is inside
    assert(polygonmask != NULL);
    while ((leafIterator != leafEnd) &&
           (!polygonmask->pointInPolygon(BoundaryXY<T>::getX(*it),
                                         BoundaryXY<T>::getY(*it)) ||
            already.end() != std::find(already.begin(), already.end(), &(*it)))
          )
    {
      ++it;
      advanceToNextLeaf();
    }
  }

  return *this;
}

template<typename T>
typename SmartQuadtree<T>::iterator::reference
SmartQuadtree<T>::iterator::operator*()
{ return it.operator*(); }

template<typename T>
typename SmartQuadtree<T>::iterator::pointer
SmartQuadtree<T>::iterator::operator->()
{ return it.operator->(); }

template<typename T> bool
SmartQuadtree<T>::iterator::operator==(
    const typename SmartQuadtree<T>::iterator& rhs) const
{
  if (it == itEnd) return leafIterator == rhs.leafIterator;
  return (it == rhs.it) && (leafIterator == rhs.leafIterator);
}

template<typename T> bool
SmartQuadtree<T>::iterator::operator!=(
    const typename SmartQuadtree<T>::iterator& rhs) const
{ return !(*this == rhs); }
/*
template<typename T>
void SmartQuadtree<T>::iterate(bool (*apply)(T&))
{
  if (children[0] != NULL)
  {
    children[0]->iterate(apply);
    children[1]->iterate(apply);
    children[2]->iterate(apply);
    children[3]->iterate(apply);
    return;
  }

  typename std::list<T>::iterator it = points.begin(), ie = points.end();
  while ( it != ie )
  {
    if (already.end() == std::find(already.begin(), already.end(), &(*it)))
      if (apply(*it))
        if (!b.contains(&(*it)))
        {
          // Computing the proper neighbour is probably slower than finding it
          // from the ancestor node...
          SmartQuadtree<T>* previous = ancestor->where[&(*it)];
          assert (previous != NULL);
          ancestor->insert(*it);
          SmartQuadtree<T>* current = ancestor->where[&(*it)];
          assert (current != NULL);
          if (current->location > previous->location)
            current->already.push_back(&(*it));
          it = points.erase(it);
        }
        else ++it;
      else ++it;
    else ++it;
  }
  already.clear();
}
*/
template<typename T>
void SmartQuadtree<T>::iterate(void (*apply)(T&, T&))
{
  if (children[0] != NULL)
  {
    children[0]->iterate(apply);
    children[1]->iterate(apply);
    children[2]->iterate(apply);
    children[3]->iterate(apply);
  }

  std::vector<const T*> neighbours;

  SmartQuadtree* nb;
  for (size_t i = 0; i < 4; ++i)
    if (delta[i] < 1) {
      nb = samelevel(i);
      assert (nb != NULL);
      typename list<T>::const_iterator it = nb->getPoints().begin();
      for ( ; it != nb->getPoints().end(); ++it)
        neighbours.push_back(&(*it));
    }
  for (size_t i = 4; i < 8; ++i)
    if (delta[i] < 0)
    {
      nb = samelevel(i);
      assert (nb != NULL);
      typename list<T>::const_iterator it = nb->getPoints().begin();
      for ( ; it != nb->getPoints().end(); ++it)
        neighbours.push_back(&(*it));
    }

  typename list<T>::iterator it = points.begin(), ie = points.end();
  for ( ; it != ie ; ++it)
  {
    typename list<T>::iterator jt = it; ++jt;
    for ( ; jt != ie; ++jt)
      apply(*it, *jt);
    for (int j = 0; j < neighbours.size(); ++j)
      apply(*it, const_cast<T&>(*neighbours[j]));
//       apply(*it, *neighbours[j]);
  }

}

template<typename T>
void SmartQuadtree<T>::iterate(const PolygonMask& m, void (*apply)(T&, T&))
{
  PolygonMask clip = m.clip(b);
  if (clip.getSize() < 3) return;

  if (children[0] != NULL)
  {
    children[0]->iterate(m, apply);
    children[1]->iterate(m, apply);
    children[2]->iterate(m, apply);
    children[3]->iterate(m, apply);
  }

  std::vector<const T*> neighbours;

  SmartQuadtree* nb;
  PolygonMask clip_nb(clip);
  for (size_t i = 0; i < 4; ++i)
    if (delta[i] < 1) {
      nb = samelevel(i);
      clip_nb = m.clip(nb->b);
      if (clip_nb.getSize() < 3) continue;
      typename list<T>::const_iterator it = nb->getPoints().begin(),
        itend = nb->getPoints().end();
     if (b.coveredByPolygon(clip) == 4)
        for (; it!=itend; ++it)
          neighbours.push_back(&(*it));
      else
        for (; it!=itend; ++it)
          if (m.pointInPolygon(BoundaryXY<T>::getX(*it),
                               BoundaryXY<T>::getY(*it)))
            neighbours.push_back(&(*it));
    }
  for (size_t i = 4; i < 8; ++i)
    if (delta[i] < 0)
    {
      nb = samelevel(i);
      clip_nb = m.clip(nb->b);
      if (clip_nb.getSize() < 3) continue;
      typename list<T>::const_iterator it = nb->getPoints().begin(),
        itend = nb->getPoints().end();
      if (b.coveredByPolygon(clip) == 4)
        for (; it!=itend; ++it)
            neighbours.push_back(&(*it));
      else
        for (; it!=itend; ++it)
          if (m.pointInPolygon(BoundaryXY<T>::getX(*it),
                               BoundaryXY<T>::getY(*it)))
            neighbours.push_back(&(*it));
    }

  typename list<T>::iterator it = points.begin(), ie = points.end();
  /*
   * Here we can do things fast
   */
  if (b.coveredByPolygon(clip) == 4)
    for ( ; it != ie ; ++it)
    {
      typename list<T>::iterator jt = it; ++jt;
      for ( ; jt != ie; ++jt)
        apply(*it, *jt);
      for (int j = 0; j < neighbours.size(); ++j)
        apply(*it, const_cast<T&>(*neighbours[j]));
//         apply(*it, *neighbours[j]);
    }
  else
    /*
     * Or we have to check more...
     */
    for ( ; it != ie ; ++it)
    {
      typename list<T>::iterator jt = it; ++jt;
      if (m.pointInPolygon(BoundaryXY<T>::getX(*it),
                           BoundaryXY<T>::getY(*it)))
      {
        for ( ; jt != ie; ++jt)
          if (m.pointInPolygon(BoundaryXY<T>::getX(*it),
                               BoundaryXY<T>::getY(*it)))
            apply(*it, *jt);

      for (int j = 0; j < neighbours.size(); ++j)
        apply(*it, const_cast<T&>(*neighbours[j]));
//         apply(*it, *neighbours[j]);
    }
  }

}
/*
template<typename T>
void SmartQuadtree<T>::iteratebyn(void (*apply)(T&, T&),
                                  void (*applybyn)(T&, T&),
                                  unsigned char n)
{
  if (children[0] != NULL)
  {
    children[0]->iteratebyn(apply, applybyn, n);
    children[1]->iteratebyn(apply, applybyn, n);
    children[2]->iteratebyn(apply, applybyn, n);
    children[3]->iteratebyn(apply, applybyn, n);
  }

  std::vector<T*> neighbours;
  typename list<T>::const_iterator it = points.begin(), ie = points.end();

  SmartQuadtree* nb;
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

  for (list<void*>::iterator itL = points.begin();
       itL != points.end(); itL++)
  {
    list<void*>::iterator jtL = itL;
    jtL++;
    for ( ; jtL != points.end(); ++jtL)
      apply(*itL, *jtL);
    int j=0;
    for ( ; j < int(neighbours.size()-n+1); j+=n)
      applybyn(*itL, &(neighbours[j]));
    for( ; j < neighbours.size(); j++)
      apply(*itL, neighbours[j]);
  }
}
*/

template<typename T>
unsigned long SmartQuadtree<T>::getDataSize() const
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

template<typename T>
unsigned char SmartQuadtree<T>::getDepth() const
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

template<typename T>
typename SmartQuadtree<T>::const_iterator MaskedQuadtree<T>::begin() const
{
  return typename SmartQuadtree<T>::const_iterator(
      quadtree.leaves.begin(), quadtree.leaves.end(), polygonmask);
}

template<typename T>
typename SmartQuadtree<T>::const_iterator MaskedQuadtree<T>::end() const
{
  return typename SmartQuadtree<T>::const_iterator(
      quadtree.leaves.end(), quadtree.leaves.end(), polygonmask);
}

template<typename T>
typename SmartQuadtree<T>::iterator MaskedQuadtree<T>::begin()
{
  return typename SmartQuadtree<T>::iterator(
      quadtree.leaves.begin(), quadtree.leaves.end(), polygonmask);
}

template<typename T>
typename SmartQuadtree<T>::iterator MaskedQuadtree<T>::end()
{
  return typename SmartQuadtree<T>::iterator(
      quadtree.leaves.end(), quadtree.leaves.end(), polygonmask);
}
