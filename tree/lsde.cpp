#include "tree.h"

GRBModel Tree::lsde(GRBEnv &env)
{
  GRBModel model(env);
  unordered_map<int, vvar> lsde_vars;

  queue<int> nodes;
  nodes.push(0);

  while (not nodes.empty())
  {
    int node = nodes.front();
    for (int child : d_children[node])
      nodes.push(child);
    nodes.pop();

    vvar parent_vars;
    if (node > 0)
      parent_vars = lsde_vars[d_ancestor[node]];

    lsde_vars[node] = add_to_lsde(node, model, parent_vars);
  }

  model.update();
  return model;
}