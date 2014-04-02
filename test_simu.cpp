#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cfloat>
#include <cstdlib>

#include <iostream>
#include <vector>

#include "quadtree.h"

GLint height = 600;
GLint width = 900;
ExtendedQuadtree* q = NULL;
float center_x, center_y;
float zoom;
int mouse_x, mouse_y, mouse_b, mouse_s;


struct Point {

  float x, y, vx, vy;
  bool draw;
  std::list<Point*> neighbours;

  Point(float x, float y) : x(x), y(y), draw(false)
  {
    vx = (((float) rand())/ (float) RAND_MAX) - 0.5;
    vy = (((float) rand())/ (float) RAND_MAX) - 0.5;
  }

  void iterate() {
    x += vx;
    y += vy;
    if (x < 0) { x = -x; vx = -vx; }
    if (y < 0) { y = -y; vy = -vy; }
    if (x > width) { vx = -vx; x += vx; }
    if (y > height) { vy = -vy; y += vy; }
    draw = false;
    neighbours.clear();
  }

  float distance2 (Point p) {
    return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
  }

};

float getX(void* p) { return ((Point*) p)->x; }
float getY(void* p) { return ((Point*) p)->y; }

bool limitation(Boundary* b) {
  float sq_size = b->norm_infty();
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
    " (0x" << std::hex << e.location << ") #" << std::dec << e.level << " -> ";
  std::list<void*>::const_iterator it = e.points.begin(),
    ie = e.points.end();
  for ( ; it != ie; ++it)
    os << *((Point*)(*it)) << " ";
  os << std::endl;
  if (NULL != e.getChild(0))
    os <<
      *(e.getChild(0)) << *(e.getChild(1)) <<
      *(e.getChild(2)) << *(e.getChild(3));
  for (size_t i=0; i<e.level; ++i) os << "  ";
  os << "}" << std::endl;
  return os;
}

bool printQuadtree(void* it)
{
  Point* p = (Point*)(it);
  if (p->draw)
    glColor3ub(185, 95, 86);
  else
    glColor3ub(65, 65, 65);

  glBegin(GL_POINTS);
  {
    glVertex3f(p->x, p->y, 1);
  }
  glEnd();

  std::list<Point*>::iterator pit = p->neighbours.begin(),
    pend = p->neighbours.end();
  for ( ; pit != pend; ++pit)
  {
    glColor3ub(185, 95, 86);
    glBegin(GL_LINES);
    {
      glVertex3f(p->x, p->y, 1);
      glVertex3f((*pit)->x, (*pit)->y, 1);
    }
    glEnd();
  }
  // does not modify the quadtree
  return false;
}

bool movePoints(void* p) {
  ((Point*) p)->iterate();
  // may modify the quadtree
  return true;
}

void conflict(void* pv1, void* pv2)
{
  Point *p1 = (Point*) pv1, *p2 = (Point*) pv2;
  if (p1->distance2(*p2) < 16.)
  {
    p1->draw = true;
    p2->draw = true;
    p1->neighbours.push_back(p2);
  }
}


void onDisplay(void)
{
  glClearColor(1., 1., 1., 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  glViewport(0,0,width,height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0,(GLfloat)width/(GLfloat)height,50.0f,15000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(center_x, center_y, zoom, // position
            center_x, center_y, 0., // view
            0.0, 1.0, 0.0); // upvector

  glPointSize(4.f);
  glLineWidth(1.f);
  q->iterate(printQuadtree);

  // The lines for the first subdivision of the quadtree
  glColor3ub(200, 200, 200);
  glLineWidth(.000000000001f);
  glBegin(GL_LINES);
  {
    glVertex3f(0, height/2, 0);
    glVertex3f(width, height/2, 0);
  }
  glEnd();
  glBegin(GL_LINES);
  {
    glVertex3f(width/2, 0, 0);
    glVertex3f(width/2, height, 0);
  }
  glEnd();


  glutSwapBuffers();

}

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

void onKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
  case 'k': // up
    center_y += 30.;
    break;
  case 'j': // down
    center_y -= 30.;
    break;
  case 'l': // right
    center_x += 30.;
    break;
  case 'h': // left
    center_x -= 30.;
    break;
  case '-': // zoom out;
    zoom *= 1.1;
    break;
  case '+': // zoom in;
    zoom *= 0.9;
    if (zoom < 100) zoom = 100;
    break;
  case 'q':
  case 27: // ESC
    std::cout << std::endl;
    exit (EXIT_SUCCESS);
  }

}
unsigned long frameCount, fps;
unsigned long currentTime, previousTime;

void calculateFPS()
{
    frameCount++;
    //  Get the number of milliseconds since glutInit called
    currentTime = glutGet(GLUT_ELAPSED_TIME);
    int timeInterval = currentTime - previousTime;

    if (timeInterval > 1000)
    {
        fps = frameCount / (timeInterval / 1000.0f);
        previousTime = currentTime;
        frameCount = 0;
    }
    std::cout << "\rfps: " << fps <<
      " size: " << q->getDataSize() <<
      " depth: " << q->getDepth() << std::flush;
}

void onIdle(void) {

  q->iterate(movePoints);
  q->iterate(conflict);

  calculateFPS();
  glutPostRedisplay();
}


int main(int argc, char* argv[])
{
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(width, height);
  glutInitWindowPosition(0,0);
  glutCreateWindow("Distance computation with quadtrees");

  glEnable (GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);

  center_x = width/2;
  center_y = height/2;
  zoom = 600;

  glutKeyboardFunc(onKeyboard);
  glutMouseFunc(onClick);
  glutMotionFunc(onMotion);
  glutDisplayFunc(onDisplay);
  glutIdleFunc(onIdle);

  q = new ExtendedQuadtree((float) width/2., (float) height/2.,
                           (float) width/2., (float) height/2., 16);

  q->setXYFcts(getX, getY);
  q->setLimitation(limitation);

  for (int i = 0; i< 10000; ++i)
    q->insert(new Point((((float) rand())/ (float) RAND_MAX) * width,
                        (((float) rand())/ (float) RAND_MAX) * height));


  glutMainLoop();
  return EXIT_SUCCESS;

}
