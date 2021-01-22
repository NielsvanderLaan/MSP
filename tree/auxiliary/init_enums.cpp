#include "../tree.h"

void Tree::init_enums(int node, Solution const &sol)
{
  for (int child : d_children[node])
  {
    Master &master = d_masters[child];
    master.update(sol);
    master.solve_mip();

    Solution copy = sol;
    copy.extend(master.mip_xvals(), master.mip_theta());

    d_enumerators[child].set_mp(sol);
    d_enumerators[child].add_point(copy);
  }
}