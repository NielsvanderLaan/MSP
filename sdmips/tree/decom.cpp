#include "tree.h"

void Tree::decom(GRBEnv &env)
{
  assert(d_masters.size() == 0);
  d_masters.reserve(d_nodes.size());
  d_enumerators.reserve(d_nodes.size());
  d_fenchel.reserve(d_nodes.size());

  for (size_t node = 0; node != d_nodes.size(); ++node)
  {
    bool leaf = is_leaf(node);
    vector<int> path = path_to(node);

    d_masters.emplace_back(Master{d_nodes[node], leaf,env});
    d_enumerators.emplace_back(Enumerator{d_nodes, path,path.size() - 1, 0, leaf,env});
    d_fenchel.emplace_back(Enumerator {d_nodes, path, path.size(), 0, leaf, env});
  }

}