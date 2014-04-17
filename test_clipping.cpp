#include "quadtree.h"
#include "logger.h"

class Test_PolygonMask {
public:
  static void RunTest_PolygonMask(Logger& log) ;
};

void Test_PolygonMask::RunTest_PolygonMask(Logger& log)
{

  Boundary b(0., 0., 10., 10.);

  {
    log.message(__LINE__, "Test of triangle intersecting with bbox");
    std::vector<float> polyX, polyY;

    polyX.push_back(-5.);   polyY.push_back(-20.);
    polyX.push_back(-15.);  polyY.push_back(5.);
    polyX.push_back(5.);    polyY.push_back(5.);

    PolygonMask m(polyX, polyY, 3);
    PolygonMask clip = m.clip(b);

    log.testint(__LINE__, clip.polyX.size(), 5, "clip.size()");

    log.testint(__LINE__, clip.polyX[0], -1.,  "clip.polyX[0]");
    log.testint(__LINE__, clip.polyX[1], -9.,  "clip.polyX[1]");
    log.testint(__LINE__, clip.polyX[2], -10., "clip.polyX[2]");
    log.testint(__LINE__, clip.polyX[3], -10., "clip.polyX[3]");
    log.testint(__LINE__, clip.polyX[4], 5.,   "clip.polyX[4]");
    log.testint(__LINE__, clip.polyY[0], -10., "clip.polyY[0]");
    log.testint(__LINE__, clip.polyY[1], -10., "clip.polyY[1]");
    log.testint(__LINE__, clip.polyY[2], -7.,  "clip.polyY[2]");
    log.testint(__LINE__, clip.polyY[3], 5.,   "clip.polyY[3]");
    log.testint(__LINE__, clip.polyY[4], 5.,   "clip.polyY[4]");
  }

  {
    log.message(__LINE__, "Test of polygon surrounding bbox");
    std::vector<float> polyX, polyY;

    polyX.push_back(-15.);  polyY.push_back(-15.);
    polyX.push_back(-15.);  polyY.push_back(15.);
    polyX.push_back(15.);   polyY.push_back(15.);
    polyX.push_back(15.);   polyY.push_back(-15.);

    PolygonMask m(polyX, polyY, 4);
    PolygonMask clip = m.clip(b);

    log.testint(__LINE__, clip.polyX.size(), 4, "clip.size()");

    log.testint(__LINE__, clip.polyX[0], 10.,  "clip.polyX[0]");
    log.testint(__LINE__, clip.polyX[1], 10.,  "clip.polyX[1]");
    log.testint(__LINE__, clip.polyX[2], -10., "clip.polyX[2]");
    log.testint(__LINE__, clip.polyX[3], -10., "clip.polyX[3]");
    log.testint(__LINE__, clip.polyY[0], 10.,  "clip.polyY[0]");
    log.testint(__LINE__, clip.polyY[1], -10., "clip.polyY[1]");
    log.testint(__LINE__, clip.polyY[2], -10., "clip.polyY[2]");
    log.testint(__LINE__, clip.polyY[3], 10.,  "clip.polyY[3]");
  }

  {
    log.message(__LINE__, "Test of polygon sharing an edge with bbox");
    std::vector<float> polyX, polyY;

    polyX.push_back(-10.);  polyY.push_back(-5.);
    polyX.push_back(-10.);  polyY.push_back(5.);
    polyX.push_back(15.);   polyY.push_back(15.);
    polyX.push_back(15.);   polyY.push_back(-15.);

    PolygonMask m(polyX, polyY, 4);
    PolygonMask clip = m.clip(b);

    log.testint(__LINE__, clip.polyX.size(), 6, "clip.size()");

    log.testint(__LINE__, clip.polyX[0], 10.,  "clip.polyX[0]");
    log.testint(__LINE__, clip.polyX[1], 10.,  "clip.polyX[1]");
    log.testint(__LINE__, clip.polyX[2], 2.,   "clip.polyX[2]");
    log.testint(__LINE__, clip.polyX[3], -10., "clip.polyX[3]");
    log.testint(__LINE__, clip.polyX[4], -10., "clip.polyX[4]");
    log.testint(__LINE__, clip.polyX[5], 2.,   "clip.polyX[5]");
    log.testint(__LINE__, clip.polyY[0], 10.,  "clip.polyY[0]");
    log.testint(__LINE__, clip.polyY[1], -10., "clip.polyY[1]");
    log.testint(__LINE__, clip.polyY[2], -10., "clip.polyY[2]");
    log.testint(__LINE__, clip.polyY[3], -5.,  "clip.polyY[3]");
    log.testint(__LINE__, clip.polyY[4], 5.,   "clip.polyY[4]");
    log.testint(__LINE__, clip.polyY[5], 10.,  "clip.polyY[5]");
  }

  {

    log.message(__LINE__, "One point on an edge, the next one out");
    Boundary subbox(225., 225., 225., 225.);
    std::vector<float> polyX, polyY;

    polyX.push_back(225.);  polyY.push_back(150.);
    polyX.push_back(225.);  polyY.push_back(300.);
    polyX.push_back(450.);   polyY.push_back(450.);
    polyX.push_back(675.);   polyY.push_back(450.);
    polyX.push_back(450.);   polyY.push_back(150.);

    PolygonMask m(polyX, polyY, 4);
    PolygonMask clip = m.clip(subbox);

    log.testint(__LINE__, clip.polyX.size(), 4, "clip.size()");

    log.testint(__LINE__, clip.polyX[0], 225.,  "clip.polyX[0]");
    log.testint(__LINE__, clip.polyX[1], 225.,  "clip.polyX[1]");
    log.testint(__LINE__, clip.polyX[2], 450.,   "clip.polyX[2]");
    log.testint(__LINE__, clip.polyX[3], 450., "clip.polyX[3]");
    log.testint(__LINE__, clip.polyY[0], 150.,  "clip.polyY[0]");
    log.testint(__LINE__, clip.polyY[1], 300., "clip.polyY[1]");
    log.testint(__LINE__, clip.polyY[2], 450., "clip.polyY[2]");
    log.testint(__LINE__, clip.polyY[3], 150.,  "clip.polyY[3]");
  }

}

int main()
{
  Logger log(__FILE__);
  Test_PolygonMask::RunTest_PolygonMask(log);
  return log.reportexit();
}
