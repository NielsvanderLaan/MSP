#include "tree.h"

void Tree::decom(GRBEnv &env)
{
  assert(d_masters.size() == 0);
  d_masters.reserve(d_nodes.size());
  d_enumerators.reserve(d_nodes.size());

  for (size_t node = 0; node != d_nodes.size(); ++node)
  {
    bool leaf = is_leaf(node);
    vector<int> path = path_to(node);

    Master master{d_nodes[node], leaf,env};
    Enumerator enumerator{d_nodes, path,path.size() - 1, leaf,env};
    Enumerator fenchel{d_nodes, path, path.size(), leaf, env};

    d_masters.push_back(master);
    d_enumerators.push_back(enumerator);
    d_fenchel.push_back(fenchel);

  }
}