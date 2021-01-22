#include "../tree.h"

vector<int> Tree::path_to(int id)
{
  vector<int> ret;
  ret.reserve(d_nodes[id].d_stage);

  int node = id;
  while (node != -1)
  {
    ret.push_back(node);
    node = d_ancestor[node];
  }

  reverse(ret.begin(), ret.end());
  return ret;
}