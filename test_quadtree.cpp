#include "quadtree.h"
#include "logger.h"

#include <cfloat> // FLT_EPSILON

struct Point {
  float x, y;
  Point(float x, float y) : x(x), y(y) {}
};

float getX(void* p) { return ((Point*) p)->x; }
float getY(void* p) { return ((Point*) p)->y; }

bool limitation(Boundary* b) {
  return (b->norm_l1() < (1 + FLT_EPSILON));
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  os << "(" << p.x << ", " << p.y << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const ExtendedQuadtree& e)
{
  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "{" << std::endl;
  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "  " <<
    e.b.center_x << ", " << e.b.center_y <<
    " (0x" << std::hex << e.location << ") #" << std::dec << e.level << " [" <<
    e.ds << "," << e.de << "," << e.dn << "," << e.dw << "] -> ";
  std::list<void*>::const_iterator it = e.points.begin(),
    ie = e.points.end();
  for ( ; it != ie; ++it)
    os << *((Point*)(*it)) << " ";
  os << std::endl;
  if (NULL != e.children[0])
    os <<
      *(e.children[0]) << *(e.children[1]) <<
      *(e.children[2]) << *(e.children[3]);
  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "}" << std::endl;
  return os;
}

int main()
{
  Logger log(__FILE__);

  ExtendedQuadtree q(0., 0., 4., 4., 4);
  q.setXYFcts(getX, getY);
  q.setLimitation(limitation);

  q.insert(new Point(1.  , 1.  ));
  q.insert(new Point(1.  , 2.  ));
  q.insert(new Point(-2. , 1.  ));
  q.insert(new Point(0.  , 2.  ));
  q.insert(new Point(0.1 , 2.  ));
  q.insert(new Point(1.  , -1. ));
  q.insert(new Point(1.  , 3.  ));
  q.insert(new Point(-2. , 2.  ));
  q.insert(new Point(1.2 , 1.3 ));
  q.insert(new Point(0.1 , 0.3 ));
  q.insert(new Point(0.1 , 0.1 ));
  q.insert(new Point(0.1 , 0.2 ));

  log.message(__LINE__, "Neighbour test in quadtree");

  // Original quadtree at this point
  //
  //                          |            |
  //                          |    0xe     |     0xf
  //                          |            |
  //            0x2           |------------|-----------
  //                          | 0x32  0x33 |
  //                          |            |     0xd
  //                          | 0x30  0x31 |
  // -------------------------|------------|-----------
  //                          |
  //                          |
  //                          |
  //            0x0           |           0x1
  //                          |
  //                          |
  //                          |
  //

  ExtendedQuadtree* m = q.getQuadrant(0x30, 3);

  log.testhex(__LINE__, m->samelevel(NORTH)->getLocation(), 0x32,
              "m->samelevel(NORTH)->getLocation()");
  log.testhex(__LINE__, m->samelevel(WEST)->getLocation(), 0x02,
              "m->samelevel(WEST)->getLocation()");
  log.testhex(__LINE__, m->samelevel(SOUTHWEST)->getLocation(), 0x00,
              "m->samelevel(SOUTHWEST)->getLocation()");
  log.testhex(__LINE__, m->samelevel(SOUTHEAST)->getLocation(), 0x01,
              "m->samelevel(SOUTHEAST)->getLocation()");
  log.testhex(__LINE__, m->samelevel(EAST)->getLocation(), 0x31,
              "m->samelevel(EAST)->getLocation()");
  log.testhex(__LINE__, m->samelevel(NORTHEAST)->getLocation(), 0x33,
              "m->samelevel(NORTHEAST)->getLocation()");
  log.testhex(__LINE__, m->samelevel(NORTH)->samelevel(NORTH)->getLocation(),
              0x0e, "m->samelevel(NORTH)->samelevel(NORTH)->getLocation()");

  log.message(__LINE__, "Test of level differences");

  log.testint(__LINE__, m->ds, -2, "m->ds");
  log.testint(__LINE__, m->samelevel(SOUTH)->dn, 1, "m->samelevel(SOUTH)->dn");
  log.testint(__LINE__, m->dn, 0, "m->dn");
  log.testint(__LINE__, m->samelevel(NORTH)->ds, 0, "m->samelevel(NORTH)->ds");
  log.testint(__LINE__, m->de, 0, "m->de");
  log.testint(__LINE__, m->samelevel(EAST)->dw, 0, "m->samelevel(EAST)->dw");
  log.testint(__LINE__, m->dw, -2, "m->dw");
  log.testint(__LINE__, m->samelevel(WEST)->de, 1, "m->samelevel(WEST)->de");

  ExtendedQuadtree* me = q.getQuadrant(0x31, 3);
  log.testint(__LINE__, me->de, -1, "me->de");
  log.testint(__LINE__, me->samelevel(EAST)->dw, 1, "me->samelevel(EAST)->dw");

  log.message(__LINE__, "Test of level differences after further insertion");

  q.insert(new Point(-1.,  1.));
  q.insert(new Point(-1.2, 1.3));

  log.testint(__LINE__, m->dw, -1, "m->dw");
  log.testint(__LINE__, m->samelevel(WEST)->de, 1, "m->samelevel(WEST)->de");

  ExtendedQuadtree* mnn = q.getQuadrant(0xe, 2);
  log.testint(__LINE__, mnn->dw, 0, "m->dw");
  log.testint(__LINE__, mnn->samelevel(WEST)->de, 0, "mnn->samelevel(WEST)->de");

  q.insert(new Point(-0.7, 0.3));
  q.insert(new Point(-0.4, 0.3));
  q.insert(new Point(-0.1, 0.6));

  ExtendedQuadtree* mw = q.getQuadrant(0x25, 3);
  log.testint(__LINE__, mw->ds, -2, "mw->ds");
  log.testint(__LINE__, mw->de, 0, "mw->de");
  log.testint(__LINE__, mw->dn, 0, "mw->dn");
  log.testint(__LINE__, mw->dw, 0, "mw->dw");

  log.testint(__LINE__, m->dw, 0, "m->dw");

  return log.reportexit();
// TODO
// test d'intégration
// passage à octree
}
