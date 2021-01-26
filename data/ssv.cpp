#include "instances.h"

Tree ssv()
{
  NodeData root {1,
                 1,
                 -320,
                 vec{-1.5, -4},
                 vec{0.0, 0.0},
                 vec{5.0, 5.0},
                 sp_mat(),
                 sp_mat(),
                 vec(),
                 Col<char>{GRB_INTEGER, GRB_INTEGER},
                 Col<char>()};

  mat tm = {{2.0/3, 1.0/3}, {1.0/3, 2.0/3}};
  mat rm = {{2, 6}, {3, 1}, {4, 3}, {5, 2}};

  sp_mat test (rm);

  NodeData sub {2,
                1.0/121,
                0.0,
                vec{-16, -19, -23, -28},
                vec{0.0, 0.0, 0.0, 0.0},
                vec{GRB_INFINITY, GRB_INFINITY, GRB_INFINITY, GRB_INFINITY},
                sp_mat(rm),
                sp_mat(tm),
                vec{5, 5},
                Col<char>{GRB_INTEGER, GRB_INTEGER, GRB_INTEGER, GRB_INTEGER},
                Col<char>{GRB_LESS_EQUAL, GRB_LESS_EQUAL}};

  Tree tree;
  tree.add_node(root);

  for (size_t s1 = 0; s1 != 11; ++s1)
  {
    for (size_t s2 = 0; s2 != 11; ++s2)
    {
      sub.d_rhs = {5.0 + s1, 5.0 + s2};
      tree.add_node(sub, 0);
    }
  }

  return tree;
}
