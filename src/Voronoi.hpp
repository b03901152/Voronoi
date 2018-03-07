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
    int x_new =
        (x - y) / 2; // for even number request, NOT divide 2 when transform
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

class FPoint {
public:
  FPoint(const FPoint &p) : _x(p._x), _y(p._y) {}
  FPoint(const Point &p) : _x(p.x), _y(p.y) {}
  FPoint(float x, float y) : _x(x), _y(y) {}
  float x() const { return _x; }
  float y() const { return _y; }
  float mDist(float x, float y) { return abs(x - _x) + abs(y - _y); }
  float mDist(const Point &i) const { return abs(i.x - _x) + abs(i.y - _y); }

  friend ostream &operator<<(ostream &os, const FPoint &p) {
    os << "(" << p.x() << ", " << p.y() << ")";
    return os;
  }

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
    // cerr << "site " << point << endl;
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
  int slop() {
    assert(s.x != e.x);
    return (s.y - e.y) / (s.x - e.x);
  }
  int length() { return s.dist(e); }
  bool bHor() { return s.y == e.y; }
  bool bVer() { return s.x == e.x; }
  void check() const {
    checkPoint(s);
    checkPoint(e);
  }
  int b() { return s.y - slop() * s.x; }
  float b_2() { return float(s.y - slop() * s.x) / 2; }
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

    // site (304, 96)
    // site (288, 104)
    // site (136, 32)
    // site (200, -88)
    // site (288, 64)

    // vX = {250, 250, 180, 20, 10, 40, 10, 30, 50, 70, 30, 50, 70, 140, 40,
    // 150};  vY = {-120, 120, -180, 20, 10, 40,  90,  50,
    //      70,   20,  10,   50, 70, 160, 120, 40};
    //
    // vX = {200, 20, 250, 200};
    // vY = {100, 30, 0, 60};
    //
    for (unsigned i = 0; i < vX.size(); ++i) {
      assert(bdL <= vX[i] && vX[i] <= bdR);
      assert(bdB <= vY[i] && vY[i] <= bdT);
      Point p(vX[i], vY[i]);
      vOriginPoints.push_back(p);
      if (bL1)
        p.L1_to_Linf();
      vSites.emplace_back(p, i);
    }
    int mDist = max(bdT - bdB, bdR - bdL) * 100;
    // vSites.emplace_back(Point(bdL - mDist, bdB - mDist),
    //                    vSites.size()); // dummy CAN NOT deleted
    // vSites.emplace_back(Point(bdL - mDist, bdT + mDist),
    //                    vSites.size()); // dummy CAN NOT deleted
    vSites.emplace_back(Point(bdL - mDist, (bdT + bdB) / 2),
                        vSites.size()); // dummy CAN NOT deleted
    vSites.emplace_back(Point(bdR + mDist, bdB - mDist),
                        vSites.size()); // dummy CAN NOT deleted
    vSites.emplace_back(Point(bdR + mDist, bdT + mDist),
                        vSites.size()); // dummy CAN NOT deleted

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
        // cerr << "add point " << r.point << endl << endl;
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
        Point e = boundaryCenter(f.point, suc->point, pre->point, true);
        if (suc->site.id != -1 && pre->site.id != -1) {
          Point s = vEdges[suc->bottom].s;
          // cerr << "add edge btw " << suc->site.point << " " <<
          // pre->site.point
          //     << " start from " << s << endl;
          vEdges.emplace_back(s, e, suc->site, pre->site);
        }
        reportEdge(*suc, f, e);
        updateInact(*suc);
        reportEdge(f, *pre, e);
        updateInact(*pre);
      } else {
        auto it = Y.find(Frontier(r.site));
        endEdge(*it);
        auto suc = Y.erase(it);
        auto pre = suc;
        --pre;
        if (pre->site.id != -1 && suc->site.id != -1) {
          reportEdge(*suc, *pre, it->end);
          // vEdges.emplace_back(it->end, suc->site, pre->site);
          // pre->top = suc->bottom = vEdges.size() - 1;
          updateInact(*pre);
          updateInact(*suc);
        }
      }
    }
    transEdge();
  }

  vector<FEdge> getEdges() {
    if (1) {
      vector<FEdge> v;
      Point m((bdL + bdR) / 2, (bdB + bdR) / 2);
      for (Edge e : vEdges) {
        e.checkSlop();
        e.check();
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
          const Point &tp = vOriginPoints[e.topId()];
          const Point &bp = vOriginPoints[e.bottomId()];
          // 8cerr << "p1 " << p1 << endl;
          // 8cerr << "tp " << tp << endl;
          // 8cerr << "bp " << bp << endl;
          assert(p1.mDist(tp) == p1.mDist(bp));
          assert(p2.mDist(tp) == p2.mDist(bp));
          int ts = e.topId();
          int bs = e.bottomId();
          if (e.bVer()) {
            ts = tp.x < bp.x ? e.topId() : e.bottomId();
            bs = tp.x < bp.x ? e.bottomId() : e.topId();
            // if (!(vOriginPoints[ts].x <= float(e.xl()) / 2)) {
            //  cerr << "tp " << tp << endl;
            //  cerr << "bp " << bp << endl;
            //  cerr << "p1 " << p1 << endl;
            //  cerr << "p2 " << p2 << endl;
            //  cerr << float(e.xl()) / 2 << endl;
            //}
            // assert(vOriginPoints[ts].x <= float(e.xl()) / 2);
            // assert(vOriginPoints[bs].x >= float(e.xl()) / 2);
          } else if (e.bHor()) {
            ts = tp.y > bp.y ? e.topId() : e.bottomId();
            bs = tp.y > bp.y ? e.bottomId() : e.topId();
            // if (!(vOriginPoints[ts].y >= float(e.yl()) / 2)) {
            //  cerr << "float(e.yl()) / 2 " << float(e.yl()) / 2 << endl;
            //  cerr << "vOriginPoints[ts] " << vOriginPoints[ts] << endl;
            //  cerr << "vOriginPoints[bs] " << vOriginPoints[bs] << endl;
            //  cerr << "p1 " << p1 << endl;
            //  cerr << "p2 " << p2 << endl;
            //  cerr << vOriginPoints[ts] << endl;
            //  cerr << vOriginPoints[bs] << endl;
            //}
            // assert(vOriginPoints[ts].y >= float(e.yl()) / 2);
            // assert(vOriginPoints[bs].y <= float(e.yl()) / 2);
          } else {
            assert(e.s.y == e.slop() * e.s.x + e.b());
            assert(e.e.y == e.slop() * e.e.x + e.b());
            if (tp.y < e.slop() * tp.x + e.b_2() ||
                bp.y > e.slop() * bp.x + e.b_2()) {
              swap(ts, bs);
            }
            assert(vOriginPoints[ts].y >=
                   e.slop() * vOriginPoints[ts].x + e.b_2());
            assert(vOriginPoints[bs].y <=
                   e.slop() * vOriginPoints[bs].x + e.b_2());
          }
          v.emplace_back(p1, p2, ts, bs);
        }
        assert(bdL <= e.s.x && e.s.x <= bdR);
        assert(bdL <= e.e.x && e.e.x <= bdR);
        assert(bdB <= e.s.y && e.s.y <= bdT);
        assert(bdB <= e.s.y && e.s.y <= bdT);
      }
      return v;
    } else {
      vector<FEdge> v;
      for (Edge &e : vEdges) {
        FPoint p1(e.e), p2(e.s);
        v.emplace_back(p1, p2, 0, 0);
      }
      return v;
    }
  }

  void endEdge(const Frontier &f) {
    // cerr << "endEdge " << f.point << " at " << f.end << endl;
    // cerr << vEdges[f.top].top.point << " " << vEdges[f.top].bottom.point
    //     << "start from " << vEdges[f.top].s << endl;
    // cerr << vEdges[f.bottom].top.point << " " <<
    // vEdges[f.bottom].bottom.point
    //     << "start from " << vEdges[f.top].s << endl;
    assert(f.top >= 0 && f.bottom >= 0);
    vEdges[f.top].e = f.end;
    vEdges[f.top].check();
    vEdges[f.top].exist = true;
    vEdges[f.bottom].e = f.end;
    vEdges[f.bottom].check();
    vEdges[f.bottom].exist = true;
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
      f.end = boundaryCenter(p, q, r, false);
      f.record = true;
      f.recordPoint.x = x;
      f.recordPoint.y = p.y;

      auto emp = X.emplace(f.site, f.recordPoint);
      assert(emp.second);
    }
  }

  Point boundaryCenter(Point obj, Point p2, Point p3, bool leftest) {
    Point p1 = obj;
    assert(p2.y >= p1.y && p1.y >= p3.y);
    int max_x = max({p1.x, p2.x, p3.x});
    int max_y = p2.y;
    int min_x = min({p1.x, p2.x, p3.x});
    int min_y = p3.y;

    vector<Point> v = {p1, p2, p3};
    sort(v.begin(), v.end());
    p1 = v[0];
    p2 = v[1];
    p3 = v[2];
    int d = max(max_x - min_x, max_y - min_y);
    int x, y;

    // swap(p1, p2);
    int d12 = p1.dist(p2);
    int d23 = p2.dist(p3);
    int d31 = p3.dist(p1);

    if (d12 >= d23 && d12 >= d31)
      swap(p1, p3);

    if (d31 >= d23 && d31 >= d12)
      swap(p1, p2);
    // p1 is not on the radium

    if (max_x - min_x == d) { // x diff is more,sol move on different y
      x = (max_x + min_x) / 2;
      if (p1.x == p2.x || p1.x == p3.x || p2.x == p3.x) {
        // set p1 to lonely one
        if (p1.x == p2.x)
          swap(p1, p3);
        else if (p1.x == p3.x)
          swap(p1, p2);
        assert(p2.x == p3.x);
        if (p1.x < p2.x) {
          y = min_y + d / 2;
        } else {
          y = max_y - d / 2;
        }
      } else {
        if (p1.y == min_y)
          y = min_y + d / 2;
        else
          y = max_y - d / 2;
      }
    } else {
      y = (min_y + max_y) / 2;
      if ((p1.x == min_x) || ((p1.y == p2.y || p1.y == p3.y || p2.y == p3.y) &&
                              !leftest)) // right most solution
        x = min_x + d / 2;
      else // left most solution
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
    // cerr << "reportEdge " << tf.point << " " << bf.point << "start at " << s
    //     << endl;
    const Site &ts = tf.site;
    const Site &bs = bf.site;
    assert(bs.point < ts.point);
    Point start = s;
    assert(start.x != INTMAX && start.y != INTMAX);
    //    cerr << "start from " << start << endl;
    vEdges.emplace_back(start, ts, bs);
    bf.top = tf.bottom = vEdges.size() - 1;
  }

  void checkVoronoi() {}

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
  vector<Point> vOriginPoints;
  set<Record> X;
  set<Frontier> Y;
  vector<Edge> vEdges;
  const int bdL;
  const int bdB;
  const int bdR; // boundary right
  const int bdT;
  bool bL1;
};
} // namespace voronoi
