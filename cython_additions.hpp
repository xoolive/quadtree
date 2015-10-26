#include "quadtree.h"

double BoundaryXY_getX(long&);
double BoundaryXY_getY(long&);
double size_limit;

template<>
double BoundaryXY<long>::getX(const long& p)
{ return BoundaryXY_getX(const_cast<long&>(p)); }
template<>
double BoundaryXY<long>::getY(const long& p)
{ return BoundaryXY_getY(const_cast<long&>(p)); }

SmartQuadtree<long>::iterator masked_begin(
    SmartQuadtree<long>& q, PolygonMask* p)
{ return MaskedQuadtree<long>(q, p).begin(); }

SmartQuadtree<long>::const_iterator masked_const_begin(
    SmartQuadtree<long>& q, PolygonMask* p)
{ return MaskedQuadtree<long>(q, p).begin(); }

template<>
bool BoundaryLimit<Boundary>::limitation(const Boundary& box)
{ return (box.norm_infty() < size_limit); }
