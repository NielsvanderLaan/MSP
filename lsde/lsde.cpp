#include "lsde.h"

LSDE::LSDE(GRBEnv &env, Tree &tree)
:
d_model(env),
d_tree(tree)
{
  queue<int> nodes;
  nodes.push(0);

  while (not nodes.empty())
  {
    int id = nodes.front();
    nodes.pop();
    add(id);

    for (int child : d_tree.d_children[id])
      nodes.push(child);
  }
}