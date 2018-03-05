#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <float.h>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <vector>
namespace voronoi {

using namespace std;

#define INTMAX INT_MAX / 10
#define INTMIN INT_MIN / 10

class FPoint {
public:
  FPoint(const FPoint &p) : _x(p._x), _y(p._y) {}
  FPoint(float x, float y) : _x(x), _y(y) {}
  float x() const { return _x; }
  float y() const { return _y; }

private:
  float _x;
  float _y;
};

class FEdge {
public:
  FEdge(const FPoint &p1, const FPoint &p2, int topSite, int bottomSite)
      : _p1(p1), _p2(p2), topSite(topSite), bottomSite(bottomSite) {}

  const FPoint &p1() const { return _p1; }
  const FPoint &p2() const { return _p2; }
  const FPoint &leftPoint() const { return _p1.x() < _p2.x() ? _p1 : _p2; }
  const FPoint &rightPoint() const { return _p1.x() < _p2.x() ? _p2 : _p1; }
  float xl() const { return min(_p1.x(), _p2.x()); }
  float xh() const { return max(_p1.x(), _p2.x()); }
  float yl() const { return min(_p1.y(), _p2.y()); }
  float yh() const { return max(_p1.y(), _p2.y()); }

  bool bHor() const { return abs(_p1.y() - _p2.y()) < 0.0001; }
  bool bVer() const { return abs(_p1.x() - _p2.x()) < 0.0001; }
  float slope() const {
    assert(!bVer());
    return float(_p1.y() - _p2.y()) / (_p1.x() - _p2.x());
  }
  float b() const {
    assert(_p1.y() - _p1.x() * slope() == _p2.y() - _p2.x() * slope());
    return _p1.y() - _p1.x() * slope();
  }
  float getY(int x) const {
    assert(!bVer());
    return x * slope() + b();
  }
  int topSiteId() const { return topSite; }
  int bottomSiteId() const { return bottomSite; }

private:
  const FPoint _p1;
  const FPoint _p2;
  const int topSite;    // or left
  const int bottomSite; // or right
};

struct Point {
  Point(int x, int y) : x(x), y(y) {}
  Point() : x(INTMAX), y(INTMAX) {}
  bool operator<(const Point &j) const {
    if (y == j.y)
      return x < j.x;
    return y < j.y;
  }

  bool operator==(const Point &i) const { return i.x == x && i.y == y; }
  friend ostream &operator<<(ostream &os, const Point &p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
  }

  const Point &operator=(const Point &i) {
    x = i.x;
    y = i.y;
    return *this;
  }

  int dist(const Point &i) const { return max(abs(i.x - x), abs(i.y - y)); }
  int mxDist(const Point &i) const { return abs(i.x - x); }
  int myDist(const Point &i) const { return abs(i.y - y); }
  void toward(const Point &obj, int dist) {
    if (mxDist(obj) < dist || myDist(obj) < dist) {
      cerr << "myDist(obj) " << myDist(obj) << endl;
      cerr << "mxDist(obj) " << mxDist(obj) << endl;
      cerr << "dist " << dist << endl;
      assert(mxDist(obj) >= dist);
      assert(myDist(obj) >= dist);
    }
    x += obj.x < x ? -dist : dist;
    y += obj.y < y ? -dist : dist;
  }
  void towardY(const Point &obj, int dist) {
    assert(myDist(obj) >= dist);
    y += obj.y < y ? -dist : dist;
  }

  void towardX(const Point &obj, int dist) {
    assert(mxDist(obj) >= dist);
    x += obj.x < x ? -dist : dist;
  }

  static int dist(const Point &i, const Point &j) {
    return max(abs(i.x - j.x), abs(i.y - j.y));
  }
  static Point middle(const Point &i, const Point &j) {
    return Point((i.x + j.x) / 2, (i.y + j.y) / 2);
  }

  void Linf_to_L1() {
    assert((x - y) % 2 == 0);
    assert((x + y) % 2 == 0);
    int x_new = (x - y) / 2; // for even number request, divide 2 when transform
                             // to FPoint for odd number
    int y_new = (x + y) / 2;
    x = x_new;
    y = y_new;
  }

  void L1_to_Linf() {
    int u = (x + y) * 2; // for even number request
    int v = (y - x) * 2; // for even number request
    x = u;
    y = v;
  }

  int x;
  int y;
};

struct Line {
  Line(Point p1, Point p2) : p1(p1), p2(p2) {}
  bool bVer() { return p1.x == p2.x; }
  bool bHor() { return p1.y == p2.y; }
  int length() { return p1.dist(p2); }
  int mxDist() { return p1.mxDist(p2); }
  int myDist() { return p1.myDist(p2); }
  int midX() { return (p1.x + p2.x) / 2; }
  int midY() { return (p1.y + p2.y) / 2; }
  int mLong() { return max(mxDist(), myDist()); }
  int mShort() { return min(mxDist(), myDist()); }
  Point &leftPoint() { return p1.x < p2.x ? p1 : p2; }
  Point &rightPoint() { return p1.x < p2.x ? p2 : p1; }
  Point &topPoint() { return p1.y < p2.y ? p2 : p1; }
  Point &bottomPoint() { return p1.y < p2.y ? p1 : p2; }
  Point &nearPoint(const Point &p) { return p.dist(p1) < p.dist(p2) ? p1 : p2; }
  bool dir() { return rightPoint().y > leftPoint().y; }
  Point middle() { return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2); }
  friend ostream &operator<<(ostream &os, const Line &l) {
    os << l.p1 << " " << l.p2 << endl;
    return os;
  }
  void bound(Point &p) {
    int max_x = max(p1.x, p2.x);
    int min_x = min(p1.x, p2.x);
    int max_y = max(p1.y, p2.y);
    int min_y = min(p1.y, p2.y);

    p.x = max(p.x, min_x);
    p.x = min(p.x, max_x);

    p.y = max(p.y, min_y);
    p.y = max(p.y, max_y);
  }

  Point p1;
  Point p2;
};

struct Site {
  Site(Point point, int id) : point(point), id(id) {
    assert(point.x % 2 == 0);
    assert(point.y % 2 == 0);
  }

  Point point;
  const int id;
};

struct Record {
  Record(const Site &s) : site(s), point(s.point), bAct(true) {}
  Record(const Site &s, Point point) : site(s), point(point), bAct(false) {}

  bool operator<(const Record &i) const {
    if (point == i.point)
      return bAct < i.bAct;
    if (point.x == i.point.x)
      return point.y < i.point.y;
    return point.x < i.point.x;
  }

  const Site site;
  const Point point;
  bool bAct; // is active
};

struct Edge {
  Edge(const Point &s, const Site &top, const Site &bottom)
      : s(s), top(top), bottom(bottom), exist(false) {
    checkPoint(s);
  }

  Edge(const Point &s, const Point &e, const Site &top, const Site &bottom)
      : s(s), e(e), top(top), bottom(bottom), exist(true) {
    check();
  }

  void checkPoint(const Point &p) const {
    int dt = top.point.dist(p) + top.point.dist(p);
    int db = bottom.point.dist(p) + bottom.point.dist(p);
    assert(dt == db);
  }
  void checkSlop() {
    assert(s.x == e.x || s.y == e.y || s.mxDist(e) == s.myDist(e));
  }
  int length() { return s.dist(e); }
  bool bHor() { return s.y == e.y; }
  bool bVer() { return s.x == e.x; }
  void check() const {
    checkPoint(s);
    checkPoint(e);
  }
  int topId() const { return top.id; }
  int bottomId() const { return bottom.id; }
  int xl() const { return min(s.x, e.x); }
  int xh() const { return max(s.x, e.x); }
  int yl() const { return min(s.y, e.y); }
  int yh() const { return max(s.y, e.y); }

  friend ostream &operator<<(ostream &ofs, const Edge &e) {
    ofs << e.s << " to " << e.e << " with site " << e.top.id << " site "
        << e.bottom.id;
    return ofs;
  }
  Point s;           // start
  Point e;           // end
  const Site top;    // top site id, or left if edge is vetical
  const Site bottom; // bottom site id, or right if edge is vetical
  bool exist;
};

struct Frontier {
  Frontier(const Site &s)
      : site(s), point(s.point), top(-1), bottom(-1), record(false) {}
  bool operator<(const Frontier &i) const { return point < i.point; }

  Site site;
  Point point;
  mutable int top;     // top edge id
  mutable int bottom;  // bottom edge id
  mutable bool record; // has record
  mutable Point recordPoint;
  mutable Point end;
};

class Voronoi {
public:
  Voronoi(vector<int> vX, vector<int> vY, int bdL, int bdB, int bdR, int bdT,
          bool bL1 = true)
      : sX(INTMIN), bdL(bdL * 2), bdB(bdB * 2), bdR(bdR * 2), bdT(bdT * 2),
        bL1(bL1) {
    vSites.reserve(vX.size());

    for (unsigned i = 0; i < vX.size(); ++i) {
      assert(bdL <= vX[i] && vX[i] <= bdR);
      assert(bdB <= vY[i] && vY[i] <= bdT);
      Point p(vX[i], vY[i]);
      if (bL1)
        p.L1_to_Linf();
      vSites.emplace_back(p, i);
    }
    int mDist = max(bdT - bdB, bdR - bdL) * 100;
    vSites.emplace_back(Point(bdL - mDist, bdB - mDist),
                        vSites.size()); // dummy
    vSites.emplace_back(Point(bdL - mDist, bdT + mDist),
                        vSites.size()); // dummy
    vSites.emplace_back(Point(bdR + mDist, bdB - mDist),
                        vSites.size()); // dummy
    vSites.emplace_back(Point(bdR + mDist, bdT + mDist),
                        vSites.size()); // dummy

    Y.emplace(Site(Point(INTMIN, INTMIN), -1));
    Y.emplace(Site(Point(INTMIN, INTMAX), -1));
    assert(Y.size() == 2);
    for (Site &s : vSites) {
      X.emplace(s);
    }

    while (!X.empty()) {
      Record r = *(X.begin());
      sX = r.point.x;
      X.erase(X.begin());
      if (r.bAct) {
        // cerr << "add r.site " << r.site.point << endl;
        auto it = Y.emplace(r.site);
        assert(it.second);

        auto suc = it.first;
        assert(it.first != Y.begin());
        assert(it.first != Y.end());
        suc++;
        assert(suc != Y.end());
        auto pre = it.first;
        pre--;
        assert(pre != Y.end());
        const Frontier &f = *(it.first);
        assert(suc->bottom == pre->top);
        Point e = boundaryCenter(suc->point, pre->point, f.point);
        if (suc->site.id != -1 && pre->site.id != -1) {
          Point s = vEdges[suc->bottom].s;
          vEdges.emplace_back(s, e, suc->site, pre->site);
        }
        reportEdge(*suc, f, e);
        updateInact(*suc);
        reportEdge(f, *pre, e);
        updateInact(*pre);
      } else {
        // cerr << "inact\n";
        auto it = Y.find(Frontier(r.site));
        endEdge(*it);
        auto pre2 = it;
        --pre2;
        auto suc = Y.erase(it);
        auto pre = suc;
        --pre;
        assert(pre == pre2);
        if (pre->site.id != -1 && suc->site.id != -1) {
          // cerr << "inact report\n";
          // cerr << "reportEdge " << suc->site.point << " " << pre->site.point
          //          << endl;
          vEdges.emplace_back(it->end, suc->site, pre->site);
          // cerr << "start from " << it->end << endl;
          pre->top = suc->bottom = vEdges.size() - 1;
          updateInact(*pre);
          updateInact(*suc);
        }
      }
    }
    transEdge();
  }

  vector<FEdge> getEdges() {
    if (bL1) {
      vector<FEdge> v;
      Point m((bdL + bdR) / 2, (bdB + bdR) / 2);
      for (Edge e : vEdges) {
        e.checkSlop();
        e.s.Linf_to_L1();
        e.e.Linf_to_L1();
        Line l(e.s, e.e);
        int ds = outDist(e.s);
        int de = outDist(e.e);
        if (ds > 0 && de > 0)
          continue;
        if (!e.length())
          continue;
        if (e.bHor()) {
          e.s.towardX(m, ds);
          e.e.towardX(m, de);
        } else if (e.bVer()) {
          e.s.towardY(m, ds);
          e.e.towardY(m, de);
        } else {
          if (e.s.y > bdT) {
            int d = e.s.y - bdT;
            e.s.y -= d;
            e.s.x -= l.dir() ? d : -d;
          }
          if (e.s.y < bdB) {
            int d = bdB - e.s.y;
            e.s.y += d;
            e.s.x += l.dir() ? d : -d;
          }
          if (e.s.x > bdR) {
            int d = e.s.x - bdR;
            e.s.x -= d;
            e.s.y -= l.dir() ? d : -d;
          }
          if (e.s.x < bdL) {
            int d = bdL - e.s.x;
            e.s.x += d;
            e.s.y += l.dir() ? d : -d;
          }

          if (e.e.y > bdT) {
            int d = e.e.y - bdT;
            e.e.y -= d;
            e.e.x -= l.dir() ? d : -d;
          }
          if (e.e.y < bdB) {
            int d = bdB - e.e.y;
            e.e.y += d;
            e.e.x += l.dir() ? d : -d;
          }
          if (e.e.x > bdR) {
            int d = e.e.x - bdR;
            e.e.x -= d;
            e.e.y -= l.dir() ? d : -d;
          }
          if (e.e.x < bdL) {
            int d = bdL - e.e.x;
            e.e.x += d;
            e.e.y += l.dir() ? d : -d;
          }
        }
        if (e.length()) {
          FPoint p1(float(e.s.x) / 2, float(e.s.y) / 2);
          FPoint p2(float(e.e.x) / 2, float(e.e.y) / 2);
          v.emplace_back(p1, p2, e.top.id, e.bottom.id);
        }
        // cerr << "bdL " << bdL << endl;
        // cerr << "bdR " << bdR << endl;
        // cerr << "bdB " << bdB << endl;
        // cerr << "bdT " << bdT << endl;
        assert(bdL <= e.s.x && e.s.x <= bdR);
        assert(bdL <= e.e.x && e.e.x <= bdR);
        assert(bdB <= e.s.y && e.s.y <= bdT);
        assert(bdB <= e.s.y && e.s.y <= bdT);
      }
      return v;
    }
  }

  void endEdge(const Frontier &f) {
    assert(f.top >= 0 && f.bottom >= 0);
    vEdges[f.top].e = f.end;
    vEdges[f.top].check();
    vEdges[f.top].exist = true;
    vEdges[f.bottom].e = f.end;
    vEdges[f.bottom].check();
    vEdges[f.bottom].exist = true;
    // cerr << "edge " << vEdges[f.top].top.point << " "
    //     << vEdges[f.top].bottom.point << " and edge "
    //     << vEdges[f.bottom].top.point << " " << vEdges[f.bottom].bottom.point
    //     << "end at " << f.end << endl;
  }

  Line edgeLines(const Edge &e) {
    vector<Line> v;
    Line el(e.top.point, e.bottom.point);
    int len = (el.mLong() - el.mShort()) / 2;
    Point m = el.middle();
    if (el.mxDist() > el.myDist()) {
      Point b = m, t = m;
      b.y -= len;
      t.y += len;
      v.emplace_back(b, t);
      // if (el.bHor())
      return v[0];

      Point t2 = t, b2 = b;
      if (el.dir()) {
        t2.x -= INTMAX / 2;
        t2.y += INTMAX / 2;
        b2.x += INTMAX / 2;
        b2.y -= INTMAX / 2;
      } else {
        t2.x += INTMAX / 2;
        t2.y += INTMAX / 2;
        b2.x -= INTMAX / 2;
        b2.y -= INTMAX / 2;
      }
      v.emplace_back(b2, b);
      v.emplace_back(t, t2);
      if (el.dir())
        swap(v[0], v[2]);
    } else {
      Point l = m, r = m;
      l.x -= len;
      r.x += len;
      v.emplace_back(l, r);
      // if (el.bVer())
      return v[0];
      Point l2 = l, r2 = r;
      if (el.dir()) {
        l2.x -= INTMAX / 2;
        l2.y += INTMAX / 2;
        r2.x += INTMAX / 2;
        r2.y -= INTMAX / 2;
      } else {
        l2.x -= INTMAX / 2;
        l2.y -= INTMAX / 2;
        r2.x += INTMAX / 2;
        r2.y += INTMAX / 2;
      }
      v.emplace_back(l2, l);
      v.emplace_back(r, r2);
    }
    return v[0];
  }

  void transEdge() {
    vector<Edge> tmp = vEdges;
    vEdges.clear();
    for (Edge &e : tmp) {
      if (!e.exist || e.top.id == -1 || e.bottom.id == -1 || e.length() == 0)
        continue;
      // cerr << "tran edge :" << e << endl;

      e.check();

      Point sp;
      Point ep;

      Line lSite(e.top.point, e.bottom.point);
      Point m = lSite.middle();
      int len = (lSite.mLong() - lSite.mShort()) / 2;
      int sd, ed;
      Line l(e.e, e.s);
      if (lSite.mxDist() <= lSite.myDist()) { // hor or diagonal
        if (e.s.x > e.e.x)
          swap(e.s, e.e);
        sp.x = m.x - len;
        sp.y = ep.y = m.y;
        ep.x = m.x + len;

        sd = e.s.myDist(m);
        ed = e.e.myDist(m);
      } else { // ver
        if (e.s.y > e.e.y)
          swap(e.s, e.e);
        sp.x = ep.x = m.x;
        sp.y = m.y - len;
        ep.y = m.y + len;

        sd = e.s.mxDist(m);
        ed = e.e.mxDist(m);
      }

      sp = e.s;
      ep = e.e;
      sp.toward(m, sd);
      ep.toward(m, ed);

      if (e.s.dist(sp) != 0) {
        // cerr << "final add edge " << vEdges.back() << endl;
        vEdges.emplace_back(e.s, sp, e.top, e.bottom);
        vEdges.back().checkSlop();
        vEdges.back().check();
      }

      if (sp.dist(ep) != 0) {

        vEdges.emplace_back(sp, ep, e.top, e.bottom);
        vEdges.back().checkSlop();
        vEdges.back().check();
        // cerr << "final add edge " << vEdges.back() << endl;
      }

      if (ep.dist(e.e) != 0) {
        vEdges.emplace_back(ep, e.e, e.top, e.bottom);
        vEdges.back().checkSlop();
        vEdges.back().check();
        // cerr << "final add edge " << vEdges.back() << endl;
      }
    }
  }

  void updateInact(const Frontier &f) {
    if (f.site.id == -1 || f.bottom == -1 || f.top == -1)
      return;

    // cerr << "updateInact " << f.point << endl;
    const Site &pre = vEdges[f.bottom].bottom;
    const Site &suc = vEdges[f.top].top;

    const Point &r = pre.point;
    const Point &p = f.site.point;
    const Point &q = suc.point;

    assert(q.y >= p.y && p.y >= r.y);

    if (r.x >= p.x && q.x >= p.x) {
      if (f.record) {
        auto record = X.find(Record(f.site, f.recordPoint));
        assert(record != X.end());
        X.erase(record);
      }

      int x = p.x + q.y - r.y;
      f.end = boundaryCenter(p, q, r);
      // cerr << "f.end " << f.end << endl;
      f.record = true;
      f.recordPoint.x = x;
      f.recordPoint.y = p.y;

      auto emp = X.emplace(f.site, f.recordPoint);
      assert(emp.second);
    }
  }

  Point boundaryCenter(Point p1, Point p2, Point p3) {
    int max_x = max({p1.x, p2.x, p3.x});
    int max_y = max({p1.y, p2.y, p3.y});
    int min_x = min({p1.x, p2.x, p3.x});
    int min_y = min({p1.y, p2.y, p3.y});

    int d = max(max_x - min_x, max_y - min_y);
    int x, y;

    int d12 = p1.dist(p2);
    int d23 = p2.dist(p3);
    int d31 = p3.dist(p1);

    if (d12 >= d23 && d12 >= d31)
      swap(p1, p3);

    if (d31 >= d23 && d31 >= d12)
      swap(p1, p2);
    // p1 is not on the radium

    if (max_x - min_x == d) {
      x = (max_x + min_x) / 2;
      if (p1.y == min_y)
        y = min_y + d / 2;
      else
        y = max_y - d / 2;
    } else {
      y = (min_y + max_y) / 2;
      if (p1.x == min_x)
        x = min_x + d / 2;
      else
        x = max_x - d / 2;
    }
    Point p(x, y);
    if (p.dist(p1) != p.dist(p2) || p.dist(p1) != p.dist(p3)) {
      cerr << "p1 " << p1 << endl;
      cerr << "p2 " << p2 << endl;
      cerr << "p3 " << p3 << endl;
      cerr << "p " << p << endl;
    }
    assert(p.dist(p1) == p.dist(p2));
    assert(p.dist(p1) == p.dist(p3));
    return p;
  }

  void reportEdge(const Frontier &tf, const Frontier &bf, Point s) {
    if (tf.site.id == -1 || bf.site.id == -1)
      return;
    const Site &ts = tf.site;
    const Site &bs = bf.site;
    assert(bs.point < ts.point);
    sX = max(ts.point.x, bs.point.x);
    Point start = s;
    assert(start.x != INTMAX && start.y != INTMAX);
    //    cerr << "start from " << start << endl;
    vEdges.emplace_back(start, ts, bs);
    bf.top = tf.bottom = vEdges.size() - 1;
  }
  void checkVoronoi() {
    int counter = 0;
    int w = 1000;
    int n = 1000;
    srand(0);
    while (true) {
      cerr << "counter " << counter++ << endl;
      vector<int> vx, vy, mx, my;
      vector<Point> sites;
      for (int i = 0; i < n; ++i) {
        vx.push_back(rand() % (w + 1));
        vy.push_back(rand() % (w + 1));
        sites.emplace_back(vx.back(), vy.back());
      }

      Voronoi vor(vx, vy, 0, 0, w, w);
      vector<vector<int>> minDist(w + 1);
      for (int i = 0; i < w + 1; ++i) {
        minDist[i].resize(w + 1, INT_MAX);
        for (int j = 0; j < w + 1; ++j) {
          Point p(i, j);
          for (int s = 0; s < n; ++s)
            minDist[i][j] = min(abs(p.x - sites[s].x) + abs(p.y - sites[s].y),
                                minDist[i][j]);
        }
      }
      struct P {
        P() : id(-1) {}
        P(int x, float y, int id, bool atUpper)
            : x(x), y(y), id(id), atUpper(atUpper) {}
        bool operator==(const P &p) const {
          return x == p.x && abs(y - p.y < 0.001) && id == p.id &&
                 atUpper == p.atUpper;
        }
        int x;
        float y;
        int id;
        bool atUpper;
      };

      struct PComp {
        bool operator()(const P &i, const P &j) const {
          if (i.id == j.id) {
            if (abs(i.y - j.y) < 0.0001)
              return i.atUpper < j.atUpper;
            return i.y < j.y;
          }
          return i.id < j.id;
        }
      };

      vector<multiset<P, PComp>> vSets(w + 1);
      for (const FEdge &e : vor.getEdges()) {
        // cerr << "edge from " << e.p1().x() << " " << e.p1().y() << " to "
        //     << e.p2().x() << " " << e.p2().y() << endl;
        for (int x = ceil(e.xl()); x <= floor(e.xh()); ++x) {
          assert(x >= 0 && x <= w);
          assert(x < vSets.size());
          if (e.bVer()) {
            vSets[x].emplace(x, e.yh(), e.topSiteId(), true);
            vSets[x].emplace(x, e.yl(), e.topSiteId(), false);

            vSets[x].emplace(x, e.yh(), e.bottomSiteId(), true);
            vSets[x].emplace(x, e.yl(), e.bottomSiteId(), false);
          } else {
            vSets[x].emplace(x, e.getY(x), e.topSiteId(), false);
            vSets[x].emplace(x, e.getY(x), e.bottomSiteId(), true);
          }
        }
      }

      for (const multiset<P, PComp> &s : vSets) {
        P pLow, pHigh;
        for (const P &p : s) {
          pHigh = p;
          if (pHigh.id == pLow.id) {
            assert(pHigh.y >= pLow.y);
            assert(pHigh.atUpper);
            assert(!pLow.atUpper);
            for (int y = ceil(pLow.y); y <= floor(pHigh.y); ++y) {
              assert(p.id < sites.size());
              if (minDist[p.x][y] !=
                  abs(p.x - sites[p.id].x) + abs(y - sites[p.id].y)) {
                cerr << "check point " << p.x << " " << y << endl;
                cerr << "site " << sites[p.id] << endl;
                cerr << "minDist[p.x][y] " << minDist[p.x][y] << endl;
                cerr << "P " << p.x << " " << p.y << " " << p.id << " "
                     << p.atUpper << endl;
              }
              assert(minDist[p.x][y] ==
                     abs(p.x - sites[p.id].x) + abs(y - sites[p.id].y));
            }
          }
        }
      }
    }
  }

  int outDist(const Point &p) const {
    int d = 0;
    d = max(d, p.x - bdR);
    d = max(d, bdL - p.x);
    d = max(d, p.y - bdT);
    d = max(d, bdB - p.y);
    return d;
  }

  int sX;
  vector<Site> vSites;
  vector<Point> vPoints; // origin position
  set<Record> X;
  set<Frontier> Y;
  vector<Edge> vEdges;
  const int bdL;
  const int bdB;
  const int bdR; // boundary right
  const int bdT;
  bool bL1;
}; // namespace voronoi
} // namespace voronoi
