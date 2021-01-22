#include "tree.h"

void Tree::init_decom(GRBEnv &env)
{
  assert(d_masters.size() == 0);
  d_masters.reserve(d_nodes.size());
  d_enumerators.reserve(d_nodes.size());

  for (size_t node = 0; node != d_nodes.size(); ++node)
    d_masters.push_back(Master(d_nodes[node],
                               d_children[node].size() == 0,
                               env));

  for (size_t node = 0; node != d_nodes.size(); ++node)
    d_enumerators.push_back(Enumerator(d_nodes,
                                       path_to(node),
                                       d_children[node].size() == 0,
                                       env));
}