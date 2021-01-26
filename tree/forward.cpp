#include "tree.h"

void Tree::forward(bool lp)
{
  queue<int> nodes;
  nodes.push(0);

  while (not nodes.empty())
  {
    int node = nodes.front();
    nodes.pop();

    solve_master(node, lp);
    // TODO
    // it could be that sol is not integer feasible in the problem associated with node.
    // In this case, we may forward the mip_solution (would this work?)
    // Do we also want to do this in the sddp?
    Solution sol = d_masters[node].lp_forward();

    for (int child : d_children[node])
    {
      d_masters[child].update(sol);
      if (not is_leaf(child))
        nodes.push(child);
    }
  }
}