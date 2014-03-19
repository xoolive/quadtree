#ifdef __APPLE__
#include <GLUT/glut.h>
#elif WIN32
#include <glut.h>
#else
#include <GL/glut.h>
#endif

#include <cstdlib>
#include <cfloat>
#include <iostream>

#include "quadtree.h"

struct Point {
  float x, y;
  bool draw;
  Point(float x, float y) : x(x), y(y), draw(false) {}
  float distance2 (Point p) {
    return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
  }
};

float getX(void* p) { return ((Point*) p)->x; }
float getY(void* p) { return ((Point*) p)->y; }

bool limitation(Boundary* b) {
  static float sq_size = b->norm_l2();
  return (sq_size < (16. + FLT_EPSILON));
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


GLint height = 800;
GLint width = 1200;
ExtendedQuadtree* q = NULL;

void printPoint(void* v, GLUquadricObj *pObj, int r, int g, int b)
{
  Point* p = (Point*) v;
  glPushMatrix();
  {
    glTranslatef(p->x, p->y, 0);
    glColor3ub(r, g, b);
    gluDisk(pObj, 0, 1, 10, 3);
  }
  glPopMatrix();
}

void printLine(void* v1, void* v2, GLUquadricObj* pObj)
{
  Point* p1 = (Point*) v1;
  Point* p2 = (Point*) v2;

  printPoint(p1, pObj, 185, 95, 86);
  printPoint(p2, pObj, 185, 95, 86);
  p1->draw = true;
  p2->draw = true;

  glColor3ub(185, 95, 86);
  glLineWidth(1);
  glBegin(GL_LINES);
  {
    glVertex3f(p1->x, p1->y, 1);
    glVertex3f(p2->x, p2->y, 1);
  }
  glEnd();
}

void printQuadtree(ExtendedQuadtree& e, GLUquadricObj* pObj)
{
  if (e.getPoints().size() > 0)
  {
    std::list<void*>::const_iterator it = e.getPoints().begin(),
      ie = e.getPoints().end();
    for ( ; it != ie; ++it) {
      if (!((Point*)(*it))->draw)
        printPoint(*it, pObj, 65, 65, 65);
      ((Point*)(*it))->draw = false;

    }
  }
  if (e.children[0] != NULL)
  {
    printQuadtree(*(e.children[0]), pObj);
    printQuadtree(*(e.children[1]), pObj);
    printQuadtree(*(e.children[2]), pObj);
    printQuadtree(*(e.children[3]), pObj);
  }
}

void printConflict(ExtendedQuadtree& q, GLUquadricObj* pObj, float distance)
{
  if (q.children[0] != NULL)
  {
    printConflict(*(q.children[0]), pObj, distance);
    printConflict(*(q.children[1]), pObj, distance);
    printConflict(*(q.children[2]), pObj, distance);
    printConflict(*(q.children[3]), pObj, distance);
  }

  std::vector<Point*> points;
  std::vector<float> _x, _y;
  std::list<void*>::const_iterator it = q.getPoints().begin(),
    ie = q.getPoints().end();
  for ( ; it != ie; ++it) { points.push_back((Point*) (*it)); }

  ExtendedQuadtree* nb;
  if (q.ds < 0)
  {
    nb = q.samelevel(SOUTH);
    it = nb->getPoints().begin(), ie = nb->getPoints().end();
    for ( ; it != ie; ++it) { points.push_back((Point*) (*it)); }
  }
  if (q.de < 1)
  {
    nb = q.samelevel(EAST);
    it = nb->getPoints().begin(), ie = nb->getPoints().end();
    for ( ; it != ie; ++it) { points.push_back((Point*) (*it)); }
  }
  if (q.dn < 1)
  {
    nb = q.samelevel(NORTH);
    it = nb->getPoints().begin(), ie = nb->getPoints().end();
    for ( ; it != ie; ++it) { points.push_back((Point*) (*it)); }
  }
  if (q.dw < 0)
  {
    nb = q.samelevel(WEST);
    it = nb->getPoints().begin(), ie = nb->getPoints().end();
    for ( ; it != ie; ++it) { points.push_back((Point*) (*it)); }
  }

  for (int i = 0; i<points.size(); ++i)
    for (int j = i + 1; j< points.size(); ++j)
      if (points[i]->distance2(*points[j]) < distance)
        printLine(points[i], points[j], pObj);
}

float center_x, center_y;
float zoom;

void onDisplay(void)
{
  glClearColor(1., 1., 1., 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glViewport(0,0,width,height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0,(GLfloat)width/(GLfloat)height,50.0f,15000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(center_x, center_y, zoom, // position
            center_x, center_y, 0., // view
            0.0, 1.0, 0.0); // upvector

  GLUquadricObj* obj = gluNewQuadric();
  printConflict(*q, obj, 16.);
  printQuadtree(*q, obj);
  gluDeleteQuadric(obj);

  glutSwapBuffers();

}

int mouse_x, mouse_y, mouse_b, mouse_s;

void onClick(int button, int state, int x, int y)
{

  if (button == GLUT_LEFT_BUTTON)
  {
    mouse_x = x;
    mouse_y = y;
    mouse_b = button;
    mouse_s = state;
    return;
  }

  if (button == 3) // scroll up
  {
    zoom *= 0.9;
    if (zoom < 100) zoom = 100;
    return;
  }
  if (button == 4) // scroll down
  {
    zoom *= 1.1;
    return;
  }
}

void onMotion(int x, int y)
{
  if (mouse_b == GLUT_LEFT_BUTTON)
  {
    center_x -= x - mouse_x;
    center_y += y - mouse_y;
    mouse_x = x;
    mouse_y = y;
  }
}

void onIdle(void) { glutPostRedisplay(); }


int main(int argc, char* argv[])
{
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(width, height);
  glutInitWindowPosition(0,0);
  glutCreateWindow("Distance simulation with quadtrees");

  glEnable (GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);


  center_x = width/2;
  center_y = height/2;
  zoom = 600;
  glutMouseFunc(onClick);
  glutMotionFunc(onMotion);
  glutDisplayFunc(onDisplay);
  glutIdleFunc(onIdle);

  q = new ExtendedQuadtree((float) width/2., (float) height/2.,
                           (float) width/2., (float) height/2., 12);

  q->setXYFcts(getX, getY);
  q->setLimitation(limitation);


  for (int i = 0; i< 4000; ++i)
    q->insert(new Point((((float) rand())/ (float) RAND_MAX) * width,
                        (((float) rand())/ (float) RAND_MAX) * height));

  glutMainLoop();
  return EXIT_SUCCESS;

}
