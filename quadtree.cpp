#include "quadtree.h"

void PolygonMask::precompute()
{
  // see http://alienryderflex.com/polygon/
  int i, j = size - 1;

  constant.resize(size);
  multiple.resize(size);

  for (i = 0; i < size; ++i)
  {
    if (polyY[j] == polyY[i])
    {
      constant[i] = polyX[i];
      multiple[i] = 0;
    }
    else
    {
      constant[i] = polyX[i] - (polyY[i] * polyX[j]) / (polyY[j] - polyY[i])
        + (polyY[i] * polyX[i]) / (polyY[j] - polyY[i]);
      multiple[i] = (polyX[j] - polyX[i])/(polyY[j] - polyY[i]);
    }

    j = i;
  }
}

PolygonMask::PolygonMask(std::vector<float> x, std::vector<float> y,
                         int size) : size(size), polyX(x), polyY(y)
{ precompute(); }

bool PolygonMask::pointInPolygon(float x, float y) const
{
  // see http://alienryderflex.com/polygon/
  int i, j = size - 1;
  bool oddNodes = false;

  for (i = 0; i < size; i++)
  {
    if ((polyY[i] < y && polyY[j] >= y) || (polyY[j] < y && polyY[i] >= y))
      oddNodes ^= (y*multiple[i] + constant[i] < x);
    j=i;
  }

  return oddNodes;
}
/*
std::ostream& operator<<(std::ostream& out, std::vector<float> x)
{
  out << "[";
  std::vector<float>::iterator it = x.begin(), ie = x.end();
  for( ; it != ie ; ++it, out<<",")
    out << *it;
  out << "]";
  return out;
}
*/

PolygonMask PolygonMask::clip(const Boundary& box) const
{
  // http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
  std::vector<float> xIn, yIn, xOut = polyX, yOut = polyY;

  std::vector<Boundary::OUTSIDE_TEST> outsideTest;
  outsideTest.push_back((Boundary::OUTSIDE_TEST) &Boundary::leftOf);
  outsideTest.push_back((Boundary::OUTSIDE_TEST) &Boundary::rightOf);
  outsideTest.push_back((Boundary::OUTSIDE_TEST) &Boundary::bottomOf);
  outsideTest.push_back((Boundary::OUTSIDE_TEST) &Boundary::upOf);

  std::vector<Boundary::INTERSECT> intersect;
  intersect.push_back((Boundary::INTERSECT) &Boundary::interLeft);
  intersect.push_back((Boundary::INTERSECT) &Boundary::interRight);
  intersect.push_back((Boundary::INTERSECT) &Boundary::interBottom);
  intersect.push_back((Boundary::INTERSECT) &Boundary::interUp);

  // for each edge of the boundary box
  for (size_t i = 0; i < 4; ++i)
  {
    xIn = xOut; yIn = yOut; xOut.clear(); yOut.clear();
    if (xIn.size() == 0) break;
    float xfrom = xIn.back(), yfrom = yIn.back();
    std::vector<float>::iterator xpoly = xIn.begin(), ypoly = yIn.begin();
    std::vector<float>::iterator xend = xIn.end();

    // for each edge of the polygon
    for ( ; xpoly != xend ; ++xpoly, ++ypoly)
    {
      if (!(box.*outsideTest[i])(*xpoly, *ypoly))
      {
        if ((box.*outsideTest[i])(xfrom, yfrom))
        {
          float x, y;
          (box.*intersect[i])(xfrom, yfrom, *xpoly, *ypoly, x, y);
          if ((x != *xpoly)||(y != *ypoly))
          {xOut.push_back(x); yOut.push_back(y);}
        }
        xOut.push_back(*xpoly); yOut.push_back(*ypoly);
      }
      else if (!(box.*outsideTest[i])(xfrom, yfrom))
      {
        float x, y;
        (box.*intersect[i])(xfrom, yfrom, *xpoly, *ypoly, x, y);
        if ((x != xfrom)||(y != yfrom))
        { xOut.push_back(x); yOut.push_back(y); }
      }
      xfrom = *xpoly; yfrom = *ypoly;
    }

  }

  return PolygonMask(xOut, yOut, xOut.size());
}


bool Boundary::contains(float x, float y)
{
  return ((x < center_x + dim_x + FLT_EPSILON) &&
          (x > center_x - dim_x - FLT_EPSILON) &&
          (y < center_y + dim_y + FLT_EPSILON) &&
          (y > center_y - dim_y - FLT_EPSILON));
}

int Boundary::coveredByPolygon(const PolygonMask& m) const
{
  int nb = 0;

  if (m.pointInPolygon(center_x + dim_x, center_y + dim_y)) ++nb;
  if (m.pointInPolygon(center_x + dim_x, center_y - dim_y)) ++nb;
  if (m.pointInPolygon(center_x - dim_x, center_y + dim_y)) ++nb;
  if (m.pointInPolygon(center_x - dim_x, center_y - dim_y)) ++nb;

  return nb;
}

void Boundary::interLeft(float x1, float y1, float x2, float y2,
                         float& xout, float& yout)
{
  xout = center_x - dim_x;
  yout = y1 + (xout - x1) / (x2 - x1) * (y2 - y1);
}

void Boundary::interRight(float x1, float y1, float x2, float y2,
                          float& xout, float& yout)
{
  xout = center_x + dim_x;
  yout = y1 + (xout - x1) / (x2 - x1) * (y2 - y1);
}

void Boundary::interBottom(float x1, float y1, float x2, float y2,
                           float& xout, float& yout)
{
  yout = center_y - dim_y;
  xout = x1 + (yout - y1) / (y2 - y1) * (x2 - x1);
}

void Boundary::interUp(float x1, float y1, float x2, float y2,
                       float& xout, float& yout)
{
  yout = center_y + dim_y;
  xout = x1 + (yout - y1) / (y2 - y1) * (x2 - x1);
}

