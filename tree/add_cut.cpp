#include "tree.h"

bool Tree::add_cut(int id, Cut &cut)
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  if (not d_masters[id].add_cut(cut))
    return false;      // no cut added

  queue<int> nodes;
  nodes.push(id);

  while (not nodes.empty())
  {
    int node = nodes.front();
    nodes.pop();
    d_enumerators[node].add_cut(cut);
    d_fenchel[node].add_cut(cut);

    for (int child : d_children[node])
      nodes.push(child);
  }

  return true;
}