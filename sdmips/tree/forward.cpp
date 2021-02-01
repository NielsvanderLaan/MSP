#include "tree.h"

void Tree::forward(bool lp)
{
  queue<int> nodes;
  nodes.push(0);

  while (not nodes.empty())
  {
    int node = nodes.front();
    nodes.pop();

    solve_master(node, true, lp);
    Solution state = d_masters[node].forward();

    for (int child : d_children[node])
    {
      d_masters[child].update(state);
      if (not is_leaf(child))
        nodes.push(child);
    }
  }
}