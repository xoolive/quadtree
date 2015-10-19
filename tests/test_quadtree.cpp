#include "quadtree.h"
#include "logger.h"

#include <cfloat> // FLT_EPSILON

struct Point {
  float x, y;
  Point(float x, float y) : x(x), y(y) {}
};

float getX(const void* p) { return ((Point*) p)->x; }
float getY(const void* p) { return ((Point*) p)->y; }

bool limitation(Boundary* b) {
  return (b->norm_l1() < (1 + FLT_EPSILON));
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  os << "(" << p.x << ", " << p.y << ")";
  return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const SmartQuadtree<Point>& e)
{
  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "{" << std::endl;

  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "  " <<
    e.b.center_x << ", " << e.b.center_y <<
    " (0x" << std::hex << e.location << ") #" << std::dec << e.level << " [" <<
    e.delta[EAST] << "," << e.delta[NORTHEAST] << "," <<
    e.delta[NORTH] << "," << e.delta[NORTHWEST] << "," <<
    e.delta[WEST] << "," << e.delta[SOUTHWEST] << "," <<
    e.delta[SOUTH] << "," << e.delta[SOUTHEAST] << "] -> ";

  std::list<Point>::const_iterator it = e.points.begin(), ie = e.points.end();
  for ( ; it != ie; ++it) os << *it << " ";
  os << std::endl;

  if (NULL != e.children[0])
    os << *(e.children[0]) << *(e.children[1]) <<
      *(e.children[2]) << *(e.children[3]);
  for (size_t i=0; i<e.level; ++i) os << "  ";

  os << "}" << std::endl;
  return os;
}

class Test_SmartQuadtree {
public:
  static void RunTest_SmartQuadtree(Logger& log) ;
};

void Test_SmartQuadtree::RunTest_SmartQuadtree(Logger& log)
{

  SmartQuadtree<Point> q(0., 0., 4., 4., 4);
  q.setXYFcts(getX, getY);
  q.setLimitation(limitation);

  q.insert(Point(1.  , 1.  ));
  q.insert(Point(1.  , 2.  ));
  q.insert(Point(-2. , 1.  ));
  q.insert(Point(0.  , 2.  ));
  q.insert(Point(0.1 , 2.  ));
  q.insert(Point(1.  , -1. ));
  q.insert(Point(1.  , 3.  ));
  q.insert(Point(-2. , 2.  ));
  q.insert(Point(1.2 , 1.3 ));
  q.insert(Point(0.1 , 0.3 ));
  q.insert(Point(0.1 , 0.1 ));
  q.insert(Point(0.1 , 0.2 ));

  log.message(__LINE__, "Tests of neighbourhood in quadtrees");

  // Original quadtree at this point
  log.message(__LINE__, "                                                   ");
  log.message(__LINE__, "                          |            |           ");
  log.message(__LINE__, "                          |    0xe     |     0xf   ");
  log.message(__LINE__, "                          |            |           ");
  log.message(__LINE__, "            0x2           |------------|-----------");
  log.message(__LINE__, "                          | 0x32  0x33 |           ");
  log.message(__LINE__, "                          |            |     0xd   ");
  log.message(__LINE__, "                          | 0x30  0x31 |           ");
  log.message(__LINE__, " -------------------------|------------|-----------");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "            0x0           |           0x1          ");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "                          |                        ");
  log.message(__LINE__, "                                                   ");

  SmartQuadtree<Point>* m = q.getQuadrant(0x30, 3);

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

  log.message(__LINE__, "");
  log.message(__LINE__, "Tests of level differences");

  // Level differences

  log.message(__LINE__, "                                                  ");
  log.message(__LINE__, " 2          2          2 | 2    2   2 | 2    2   2");
  log.message(__LINE__, "                         | -1       0 | 0        2");
  log.message(__LINE__, "                         | -1   0   0 | 0    0   2");
  log.message(__LINE__, " 2                     1 |------------|-----------");
  log.message(__LINE__, "                         | -1   0   0 | 0    0   2");
  log.message(__LINE__, "                         | -1       0 | 0        2");
  log.message(__LINE__, " 2          0          0 | -1  -1  -1 | -1  -1   2");
  log.message(__LINE__, "-------------------------|------------|-----------");
  log.message(__LINE__, " 2          0          1 | 0          1          2");
  log.message(__LINE__, "                         |                        ");
  log.message(__LINE__, "                         |                        ");
  log.message(__LINE__, " 2                     0 | 0                     2");
  log.message(__LINE__, "                         |                        ");
  log.message(__LINE__, "                         |                        ");
  log.message(__LINE__, " 2          2          2 | 2          2          2");
  log.message(__LINE__, "                                                  ");

  log.testint(__LINE__, m->delta[SOUTH], -2, "m->delta[SOUTH]");
  log.testint(__LINE__, m->samelevel(SOUTH)->delta[NORTH], 1,
              "m->samelevel(SOUTH)->delta[NORTH]");
  log.testint(__LINE__, m->delta[NORTH], 0, "m->delta[NORTH]");
  log.testint(__LINE__, m->samelevel(NORTH)->delta[SOUTH], 0,
              "m->samelevel(NORTH)->delta[SOUTH]");
  log.testint(__LINE__, m->delta[EAST], 0, "m->delta[EAST]");
  log.testint(__LINE__, m->samelevel(EAST)->delta[WEST], 0,
              "m->samelevel(EAST)->delta[WEST]");
  log.testint(__LINE__, m->delta[WEST], -2, "m->delta[WEST]");
  log.testint(__LINE__, m->samelevel(WEST)->delta[EAST], 1,
              "m->samelevel(WEST)->delta[EAST]");

  SmartQuadtree<Point>* me = q.getQuadrant(0x31, 3);
  log.testint(__LINE__, me->delta[EAST], -1, "me->delta[EAST]");
  log.testint(__LINE__, me->samelevel(EAST)->delta[WEST], 1,
              "me->samelevel(EAST)->delta[WEST]");

  log.message(__LINE__, "");
  log.message(__LINE__, "Tests of level differences in diagonal");

  log.testint(__LINE__, m->delta[SOUTHWEST], -2, "m->delta[SOUTHWEST]");
  log.testint(__LINE__, m->samelevel(SOUTHWEST)->delta[NORTHEAST], 1,
              "m->samelevel(SOUTHWEST)->delta[NORTHEAST]");

  SmartQuadtree<Point>* mne = m->samelevel(NORTHEAST);

  log.testint(__LINE__, mne->delta[NORTHEAST], -1, "mne->delta[NORTHEAST]");
  log.testint(__LINE__, mne->samelevel(NORTHEAST)->delta[SOUTHWEST], 1,
              "mne->samelevel(NORTHEAST)->delta[SOUTHWEST]");

  log.testint(__LINE__, q.getQuadrant(0, 1)->delta[NORTHEAST], 1,
              "q.getQuadrant(0, 1)->delta[NORTHEAST]");
  log.testint(__LINE__, q.getQuadrant(1, 1)->delta[NORTHWEST], 0,
              "q.getQuadrant(1, 1)->delta[NORTHWEST]");
  log.testint(__LINE__, q.getQuadrant(2, 1)->delta[SOUTHEAST], 0,
              "q.getQuadrant(2, 1)->delta[SOUTHEAST]");
  log.testint(__LINE__, q.getQuadrant(0x30, 3)->delta[SOUTHWEST], -2,
              "q.getQuadrant(3, 1)->delta[SOUTHWEST]");

  log.message(__LINE__, "");
  log.message(__LINE__, "Tests of level differences after further insertion");

  q.insert(Point(-1.,  1.));
  q.insert(Point(-1.2, 1.3));

  log.testint(__LINE__, m->delta[WEST], -1, "m->delta[WEST]");
  log.testint(__LINE__, m->samelevel(WEST)->delta[EAST], 1,
              "m->samelevel(WEST)->delta[EAST]");

  SmartQuadtree<Point>* mnn = q.getQuadrant(0xe, 2);
  log.testint(__LINE__, mnn->delta[WEST], 0, "m->delta[WEST]");
  log.testint(__LINE__, mnn->samelevel(WEST)->delta[EAST], 0,
              "mnn->samelevel(WEST)->delta[EAST]");

  log.message(__LINE__, "");
  log.message(__LINE__,
              "Tests of level differences in diagonal after further insertion");

  log.testint(__LINE__, q.getQuadrant(0x9,2)->delta[NORTHEAST], 0,
              "q.getQuadrant(0x9,2)->delta[NORTHEAST]");
  log.testint(__LINE__, q.getQuadrant(0xb,2)->delta[SOUTHEAST], 1,
              "q.getQuadrant(0xb,2)->delta[SOUTHEAST]");
  log.testint(__LINE__, q.getQuadrant(0xe,2)->delta[SOUTHWEST], 0,
              "q.getQuadrant(0xe,2)->delta[SOUTHWEST]");
  log.testint(__LINE__, q.getQuadrant(0x32,3)->delta[NORTHWEST], -1,
              "q.getQuadrant(0x32,3)->delta[NORTHWEST]");

  log.message(__LINE__, "");
  log.message(__LINE__, "Tests of level differences after further insertion");

  q.insert(Point(-0.7, 0.3));
  q.insert(Point(-0.4, 0.3));
  q.insert(Point(-0.1, 0.6));

  SmartQuadtree<Point>* mw = q.getQuadrant(0x25, 3);
  log.testint(__LINE__, mw->delta[SOUTH], -2, "mw->delta[SOUTH]");
  log.testint(__LINE__, mw->delta[EAST],   0, "mw->delta[EAST]");
  log.testint(__LINE__, mw->delta[NORTH],  0, "mw->delta[NORTH]");
  log.testint(__LINE__, mw->delta[WEST],   0, "mw->delta[WEST]");
  log.testint(__LINE__, m->delta[WEST],    0, "m->delta[WEST]");


  log.message(__LINE__, "");
  log.message(__LINE__,
              "Tests of level differences in diagonal after further insertion");

  log.testint(__LINE__, q.getQuadrant(0x27,3)->delta[NORTHEAST], -1,
              "q.getQuadrant(0x9,2)->delta[NORTHEAST]");
  log.testint(__LINE__, q.getQuadrant(0xb,2)->delta[SOUTHEAST], 1,
              "q.getQuadrant(0xb,2)->delta[SOUTHEAST]");
  log.testint(__LINE__, q.getQuadrant(0xe,2)->delta[SOUTHWEST], 1,
              "q.getQuadrant(0xe,2)->delta[SOUTHWEST]");
  log.testint(__LINE__, q.getQuadrant(0x32,3)->delta[NORTHWEST], -1,
              "q.getQuadrant(0x32,3)->delta[NORTHWEST]");
}

int main()
{
  Logger log(__FILE__);
  Test_SmartQuadtree::RunTest_SmartQuadtree(log);
  return log.reportexit();
}
