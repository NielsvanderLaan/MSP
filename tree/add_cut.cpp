#include "tree.h"

bool Tree::add_cut(int id, Cut &cut, Solution &sol)
{
  if (not d_masters[id].add_cut(cut, sol))
    return false;      // no cut added

  queue<int> nodes;
  nodes.push(id);

  while(not nodes.empty())
  {
    int node = nodes.front();
    nodes.pop();
    if (node != 0)
      d_enumerators[node].add_cut(cut);

    for (int child : d_children[node])
      nodes.push(child);
  }

  return true;
}