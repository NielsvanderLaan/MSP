#include "tree.h"

bool Tree::solve_master(int node, bool lp)
{
  Master &master = d_masters[node];

  master.solve_lp();
  if (lp)
    return true;

  while (not master.integer())
  {
    Cut cut = fenchel_cut(node);

    if (not add_cut(node, cut))
      return false;

    master.solve_lp();
  }
  return true;
}