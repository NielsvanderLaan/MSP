#include "tree.h"

void Tree::solve_master(int node, bool lp)
{
  Master &master = d_masters[node];

  master.solve_lp();
  if (lp)
    return;

  while (not master.integer())
  {
    Cut cut = fenchel_cut(node);

    if (not add_cut(node, cut))
      break;

    master.solve_lp();
  }
}