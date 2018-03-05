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
  //
  // double l = -w, b = -w, r = w, t = w;
  // vector<int> vx;
  // vector<int> vy;
  //
  // for (int i = 0; i < 10000; ++i) {
  //  vx.push_back((double)w * abs(rand()) / RAND_MAX);
  //  vy.push_back((double)w * abs(rand()) / RAND_MAX);
  //}
  // voronoi::Voronoi vor(vx, vy, l, b, r, t);
  // vor.checkVoronoi();
  // return 0;
  //
  // std::cout << "voronois done!\n";

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

void drawVoronoi() {
  using namespace std;

  double l = -w, b = -w, r = w, t = w;
  vector<int> vx;
  vector<int> vy;

  for (int i = 0; i < 100; ++i) {
    vx.push_back(int(w * ((double)rand() / RAND_MAX - 0.5)) * 2);
    vy.push_back(int(w * ((double)rand() / RAND_MAX - 0.5)) * 2);
    // cerr << vx[i] << " " << vy[i] << endl;
  }

  voronoi::Voronoi vor(vx, vy, l, b, r, t);

  glBegin(GL_QUADS);
  for (unsigned i = 0; i < vx.size(); ++i) {
    double x = vx[i];
    x /= w;
    double y = vy[i];
    y /= w;
    glVertex2f(x - 5 / w, y - 5 / w);
    glVertex2f(x + 5 / w, y - 5 / w);
    glVertex2f(x + 5 / w, y + 5 / w);
    glVertex2f(x - 5 / w, y + 5 / w);
  }
  glEnd();

  glBegin(GL_LINES);
  for (voronoi::FEdge e : vor.getEdges()) {
    e.p1().x();
    float x1 = e.p1().x();
    float y1 = e.p1().y();
    float x2 = e.p2().x();
    float y2 = e.p2().y();
    assert(x1 >= l && x1 <= r);
    assert(y1 >= l && y1 <= r);
    assert(x2 >= l && x2 <= r);
    assert(y2 >= l && y2 <= r);
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