#include "Voronoi.hpp"
#include <GL/glew.h> // Include the GLEW header file
#include <GL/glut.h> // Include the GLUT header file
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>

void display(void);
void onEF(int n);
void reshape(int width, int height);

double w = 1000;

int main(int argc, char **argv) {

  srand(time(NULL));
  using namespace std;

  // double l = -w, b = -w, r = w, t = w;
  // vector<int> vx = {100, 50};
  // vector<int> vy = {0, 0};

  // for (int i = 0; i < 100; ++i) {
  //  vx.push_back((double)w * abs(rand()) / RAND_MAX);
  //  vy.push_back((double)w * abs(rand()) / RAND_MAX);
  //}
  // voronoi::Voronoi vor(vx, vy, l, b, r, t);
  //  vor.getEdges();
  // vor.checkVoronoi();
  // return 0;

  std::cout << "voronois done!\n";

  glutInit(&argc, argv);            // Initialize GLUT
  glutInitDisplayMode(GLUT_SINGLE); // Set up a basic display buffer (only
                                    // single buffered for now)
  glutInitWindowSize(1000, 1000);   // Set the width and height of the window
  glutInitWindowPosition(0, 0);     // Set the position of the window
  glutCreateWindow(
      "You�re first OpenGL Window"); // Set the title for the window

  glutTimerFunc(100, onEF, 0);
  glutDisplayFunc(
      display); // Tell GLUT to use the method "display" for rendering

  glutReshapeFunc(
      reshape); // Tell GLUT to use the method "reshape" for reshaping

  // glutKeyboardFunc(keyPressed); // Tell GLUT to use the method "keyPressed"
  // for key presses  glutKeyboardUpFunc(keyUp); // Tell GLUT to use the
  // method "keyUp" for key up events

  glutMainLoop(); // Enter GLUT's main loop

  return 0;
}

int counter = 180;
void drawVoronoi() {
  using namespace std;
  using namespace voronoi;
  srand(counter);
  int w = 100;
  double l = 0, b = 0, r = w, t = w;

  int n = 100;
  cerr << "counter " << counter++ << endl;
  vector<int> vx, vy, mx, my;
  vector<Point> sites;
  for (int i = 0; i < n; ++i) {
    vx.push_back(rand() % (w + 1));
    vy.push_back(rand() % (w + 1));
    if (vx.back() % 2)
      ++vx.back();
    if (vy.back() % 2)
      ++vy.back();
    // cerr << "site " << vx.back() << " " << vy.back() << endl;
    sites.emplace_back(vx.back(), vy.back());
  }

  Voronoi vor(vx, vy, 0, 0, w, w);
  vector<vector<int>> minDist(w + 1);
  for (int i = 0; i < w + 1; ++i) {
    minDist[i].resize(w + 1, INT_MAX);
    for (int j = 0; j < w + 1; ++j) {
      Point p(i, j);
      for (int s = 0; s < n; ++s)
        minDist[i][j] =
            min(abs(p.x - sites[s].x) + abs(p.y - sites[s].y), minDist[i][j]);
    }
  }

  struct P {
    P(int siteId) : yh(-1), yl(1000000), siteId(siteId) {}
    bool operator==(const P &p) const {
      return abs(yh - p.yh < 0.001) && siteId == p.siteId &&
             abs(yl - p.yl < 0.001);
    }
    void setY(float y) {
      yh = max(yh, y);
      yl = min(yl, y);
    }
    float yh;
    float yl;
    const int siteId;
  };

  vector<map<int, P>> vm(w + 1);
  vector<FEdge> vf = vor.getEdges();
  for (const FEdge &e : vf) {
    // cerr << "edge " << e.p1() << " " << e.p2() << " top:" << e.topSiteId()
    //     << " bottom " << e.bottomSiteId() << endl;
    for (int x = ceil(e.xl()); x <= floor(e.xh()); ++x) {
      assert(x >= 0 && x <= w);
      assert(x < vm.size());
      if (e.bVer()) {
        auto t = vm[x].emplace(e.topSiteId(), e.topSiteId());
        t.first->second.setY(e.yh());
        t.first->second.setY(e.yl());

        auto b = vm[x].emplace(e.bottomSiteId(), e.bottomSiteId());
        b.first->second.setY(e.yh());
        b.first->second.setY(e.yl());
      } else {
        float y = e.getY(x);
        auto t = vm[x].emplace(e.topSiteId(), e.topSiteId());
        t.first->second.setY(y);

        auto b = vm[x].emplace(e.bottomSiteId(), e.bottomSiteId());
        b.first->second.setY(y);
      }
    }
  }

  for (int x = 0; x <= (int)w; ++x) {
    auto &m = vm[x];
    for (auto it : m) {
      int yBegin = ceil(it.second.yl);
      int yEnd = floor(it.second.yh);
      for (int y = yBegin; y <= yEnd; ++y) {
        cerr << "check (" << x << " " << y << endl;
        if (minDist[x][y] != abs(x - sites[it.second.siteId].x) +
                                 abs(y - sites[it.second.siteId].y)) {
          cerr << "check (" << x << " " << y << ") site ("
               << sites[it.second.siteId].x << ", " << sites[it.second.siteId].y
               << ") minDist " << minDist[x][y] << endl;
          cerr << "dist "
               << abs(x - sites[it.second.siteId].x) +
                      abs(y - sites[it.second.siteId].y)
               << endl;
        }

        assert(minDist[x][y] == abs(x - sites[it.second.siteId].x) +
                                    abs(y - sites[it.second.siteId].y));
      }
    }
  }

  glBegin(GL_QUADS);
  for (unsigned i = 0; i < vx.size(); ++i) {
    double x = vx[i];
    x /= w;
    double y = vy[i];
    y /= w;
    glVertex2f(x - 1.0 / w, y - 1.0 / w);
    glVertex2f(x + 1.0 / w, y - 1.0 / w);
    glVertex2f(x + 1.0 / w, y + 1.0 / w);
    glVertex2f(x - 1.0 / w, y + 1.0 / w);
  }
  glEnd();

  glBegin(GL_LINES);
  for (voronoi::FEdge e : vf) {
    float x1 = e.p1().x();
    float y1 = e.p1().y();
    float x2 = e.p2().x();
    float y2 = e.p2().y();
    // assert(x1 >= l && x1 <= r);
    // assert(y1 >= l && y1 <= r);
    // assert(x2 >= l && x2 <= r);
    // assert(y2 >= l && y2 <= r);
    x1 /= w;
    y1 /= w;
    x2 /= w;
    y2 /= w;
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
  }
  glEnd();
}

void display(void) {
  std::cout << "display\n";
  glLoadIdentity(); // Load the Identity Matrix to reset our drawing locations
  glTranslatef(0.0f, 0.0f, -5.0f);
  glFlush();
}

void onEF(int n) {

  glutTimerFunc(20, onEF, 0);
  glClear(GL_COLOR_BUFFER_BIT); // Clear the screen
  glClearColor(0.0f, 0.0f, 0.2f,
               1.0f); // Clear the background of our window to red

  drawVoronoi();
  glutSwapBuffers();
  std::cin.get();
  // Draw everything to the screen
}

void reshape(int width, int height) {

  glViewport(0, 0, (GLsizei)width,
             (GLsizei)height); // Set our viewport to the size of our window
  glMatrixMode(GL_PROJECTION); // Switch to the projection matrix so that we can
                               // manipulate how our scene is viewed
  glLoadIdentity(); // Reset the projection matrix to the identity matrix so
                    // that we don't get any artifacts (cleaning up)
  gluPerspective(22.5, (GLfloat)width / (GLfloat)height, 1.0,
                 100.0); // Set the Field of view angle (in degrees), the aspect
                         // ratio of our window, and the new and far planes
  glMatrixMode(GL_MODELVIEW); // Switch back to the model view matrix, so that
                              // we can start drawing shapes correctly
}