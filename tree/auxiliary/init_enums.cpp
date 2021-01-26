#include "../tree.h"

void Tree::init_enums(int node)
{
  Solution sol = d_masters[node].lp_forward();

  for (int child : d_children[node])
  {
    Master &master = d_masters[child];
    master.solve_mip();

    Solution forward = master.mip_forward();

    d_enumerators[child].set_mp(sol);
    d_enumerators[child].add_point(forward);
  }
}