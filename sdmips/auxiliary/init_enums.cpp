#include "../tree/tree.h"

void Tree::init_enums(int node)
{
  for (int child : d_children[node])
  {
    Master &sub = d_masters[child];
    d_enumerators[child].set_mp(sub.d_state);

    sub.solve_mip();
    d_enumerators[child].add_point(sub.forward(), false, true);
  }
}