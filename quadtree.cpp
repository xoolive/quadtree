#include "quadtree.h"

#include <cfloat> // DBL_EPSILON

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

PolygonMask::PolygonMask(std::vector<double> x, std::vector<double> y,
                         int size) : size(size), polyX(x), polyY(y)
{ precompute(); }

bool PolygonMask::pointInPolygon(double x, double y) const
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
std::ostream& operator<<(std::ostream& out, std::vector<double> x)
{
  out << "[";
  std::vector<double>::iterator it = x.begin(), ie = x.end();
  for( ; it != ie ; ++it, out<<",")
    out << *it;
  out << "]";
  return out;
}
*/

PolygonMask PolygonMask::clip(const Boundary& box) const
{
  // http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
  std::vector<double> xIn, yIn, xOut = polyX, yOut = polyY;

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
    double xfrom = xIn.back(), yfrom = yIn.back();
    std::vector<double>::iterator xpoly = xIn.begin(), ypoly = yIn.begin();
    std::vector<double>::iterator xend = xIn.end();

    // for each edge of the polygon
    for ( ; xpoly != xend ; ++xpoly, ++ypoly)
    {
      if (!(box.*outsideTest[i])(*xpoly, *ypoly))
      {
        if ((box.*outsideTest[i])(xfrom, yfrom))
        {
          double x, y;
          (box.*intersect[i])(xfrom, yfrom, *xpoly, *ypoly, x, y);
          if ((x != *xpoly)||(y != *ypoly))
          {xOut.push_back(x); yOut.push_back(y);}
        }
        xOut.push_back(*xpoly); yOut.push_back(*ypoly);
      }
      else if (!(box.*outsideTest[i])(xfrom, yfrom))
      {
        double x, y;
        (box.*intersect[i])(xfrom, yfrom, *xpoly, *ypoly, x, y);
        if ((x != xfrom)||(y != yfrom))
        { xOut.push_back(x); yOut.push_back(y); }
      }
      xfrom = *xpoly; yfrom = *ypoly;
    }

  }

  return PolygonMask(xOut, yOut, xOut.size());
}


bool Boundary::contains(double x, double y)
{
  return ((x < center_x + dim_x * 1.00001) &&
          (x > center_x - dim_x * 1.00001) &&
          (y < center_y + dim_y * 1.00001) &&
          (y > center_y - dim_y * 1.00001));
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

void Boundary::interLeft(double x1, double y1, double x2, double y2,
                         double& xout, double& yout)
{
  xout = center_x - dim_x;
  yout = y1 + (xout - x1) / (x2 - x1) * (y2 - y1);
}

void Boundary::interRight(double x1, double y1, double x2, double y2,
                          double& xout, double& yout)
{
  xout = center_x + dim_x;
  yout = y1 + (xout - x1) / (x2 - x1) * (y2 - y1);
}

void Boundary::interBottom(double x1, double y1, double x2, double y2,
                           double& xout, double& yout)
{
  yout = center_y - dim_y;
  xout = x1 + (yout - y1) / (y2 - y1) * (x2 - x1);
}

void Boundary::interUp(double x1, double y1, double x2, double y2,
                       double& xout, double& yout)
{
  yout = center_y + dim_y;
  xout = x1 + (yout - y1) / (y2 - y1) * (x2 - x1);
}
