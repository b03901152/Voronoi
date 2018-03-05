#include "Voronoi.hpp"
int main() {
  int n, m;
  cin >> m;
  vector<int> sweet_x, sweet_y;
  int x, y, w = 10000;

  for (int i = 0; i < m; ++i) {
    cin >> x >> y;
    sweet_x.push_back(x);
    sweet_y.push_back(y);
  }

  assert(sweet_x.size() == (unsigned)m);
  Voronoi vorornoi(sweet_x, sweet_y, 0, 0, w, w);
  assert(sweet_x.size() == (unsigned)m);
  vector<Edge> v = vorornoi.getEdges();
  assert(sweet_x.size() == (unsigned)m);

  struct Cockroach {
    Cockroach(int id, int y) : id(id), y(y), d(INTMAX) {}
    int id;
    int y;
    int d; // dist
  };

  map<int, Cockroach> mCocks;
  cin >> n;
  vector<vector<int>> output(n);
  for (int i = 0; i < n; ++i) {
    cin >> x >> y;
    mCocks.emplace(x, Cockroach(i, y));
  }
  cerr << "vor\n";

  for (const Edge &e : v) {
    cerr << "Edge " << e << endl;
    cerr << "vor\n";

    auto it = mCocks.lower_bound(e.xl());
    for (; it != mCocks.upper_bound(e.xh()); ++it) {
      cerr << "A\n";
      assert(sweet_x.size() == (unsigned)m);

      Cockroach &c = it->second;
      cerr << sweet_x.size() << endl;
      cerr << sweet_y.size() << endl;
      cerr << e.topId() << endl;
      int tx = sweet_x[e.topId()];
      int ty = sweet_y[e.topId()];
      cerr << "A\n";

      int bx = sweet_x[e.bottomId()];
      int by = sweet_y[e.bottomId()];

      cerr << "A\n";
      x = it->first;
      y = c.y;
      int td = abs(tx - x) + abs(ty - y);
      int bd = abs(bx - x) + abs(by - y);
      cerr << "B\n";
      if (td <= c.d) {
        cerr << "C\n";
        c.d = td;
        if (td < c.d)
          output[c.id].clear();
        output[c.id].push_back(e.topId());
      }
      cerr << "B\n";

      if (bd <= c.d) {
        cerr << "C\n";
        c.d = bd;
        if (bd < c.d)
          output[c.id].clear();
        output[c.id].push_back(e.bottomId());
      }
      cerr << "B\n";
    }
  };
  cerr << "output.size() " << output.size() << endl;
  for (auto &v : output) {
    for (int i : v)
      cerr << i << " @";
    cerr << endl;
  }
}