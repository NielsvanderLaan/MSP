#include "tree.h"

void Tree::forward(bool lp)
{
  queue<int> nodes;
  nodes.push(0);

  while (not nodes.empty())
  {
    int node = nodes.front();
    nodes.pop();

    //Solution state = (solve_master(node, lp) || lp) ? d_masters[node].lp_forward() : d_masters[node].mip_forward();
    solve_master(node, lp);
    Solution state = d_masters[node].lp_forward();

    for (int child : d_children[node])
    {
      d_masters[child].update(state);
      if (not is_leaf(child))
        nodes.push(child);
    }
  }
}