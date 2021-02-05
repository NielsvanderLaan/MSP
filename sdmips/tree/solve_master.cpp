#include "tree.h"

void Tree::solve_master(int node, bool affine, bool lp, bool force)
{
  Master &master = d_masters[node];

  master.solve_lp();

  if (lp)
    return;

  while (not master.integer())
  {
    Cut cut = fenchel_cut(node, affine);

    if (not add_cut(node, cut))    // mip could not be solved using cutting planes
    {
      if (force)
        master.solve_mip();           // uses Gurobi
      return;
    }

    master.solve_lp();
  }
}