#include "../tree.h"

vector<int> Tree::nvars(int id)
{
  vector<int> ret;
  vector<int> path = path_to(id);
  ret.reserve(path.size());

  for (int node : path)
    ret.push_back(d_nodes[node].nvars());

  return ret;
}
