#include "Voronoi.hpp"
int main() {
  int n, m;
  cin >> m;
  vector<int> sweet_x, sweet_y;
  int rate = 2;
  double w = 10000 * rate;
  auto trans = [&](double x) -> int { return x * rate; };
  double x, y;

  for (int i = 0; i < m; ++i) {
    cin >> x >> y;
    sweet_x.push_back(trans(x));
    sweet_y.push_back(trans(y));
  }

  assert(sweet_x.size() == (unsigned)m);
  Voronoi vorornoi(sweet_x, sweet_y, 0, 0, 2 * w, 2 * w);
  assert(sweet_x.size() == (unsigned)m);
  vector<Edge> v = vorornoi.getEdges();
  assert(sweet_x.size() == (unsigned)m);

  struct Cockroach {
    Cockroach(int id, int y) : id(id), y(y), d(INTMAX) {}
    int id;
    int y;
    int d; // dist
  };

  multimap<int, Cockroach> mCocks;
  cin >> n;
  vector<vector<int>> output(n);
  for (int i = 0; i < n; ++i) {
    cin >> x >> y;
    x = trans(x);
    y = trans(y);
    mCocks.emplace(x, Cockroach(i, y));
  }

  assert(mCocks.size() == n);
  assert(mCocks.size() == 2);
  for (const Edge &e : v) {
    auto it = mCocks.lower_bound(e.xl());
    for (; it != mCocks.upper_bound(e.xh()); ++it) {
      assert(sweet_x.size() == (unsigned)m);

      Cockroach &c = it->second;
      int tx = sweet_x[e.topId()];
      int ty = sweet_y[e.topId()];

      int bx = sweet_x[e.bottomId()];
      int by = sweet_y[e.bottomId()];

      x = it->first;
      y = c.y;
      int td = abs(tx - x) + abs(ty - y);
      int bd = abs(bx - x) + abs(by - y);

      if (td <= c.d) {
        if (td < c.d)
          output[c.id].clear();
        c.d = td;
        output[c.id].push_back(e.topId());
      }

      if (bd <= c.d) {
        if (bd < c.d)
          output[c.id].clear();
        c.d = bd;
        output[c.id].push_back(e.bottomId());
      }
    }
  };
  for (vector<int> &v : output) {
    sort(v.begin(), v.end());
    auto it = unique(v.begin(), v.end());
    v.resize(distance(v.begin(), it));
    for (int i : v)
      cerr << i + 1 << " ";
    cerr << endl;
  }
}