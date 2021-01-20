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
                 Col<char>{GRB_CONTINUOUS, GRB_CONTINUOUS},
                 Col<char>()};

  sp_mat id = eye<sp_mat>(2,2);
  mat rm = {{2, 6}, {3, 1}, {4, 3}, {5, 2}};

  sp_mat test (rm);

  NodeData sub {2,
                1.0/441,
                0.0,
                vec{-16, -19, -23, -28},
                vec{0.0, 0.0, 0.0, 0.0},
                vec{1.0, 1.0, 1.0, 1.0},
                sp_mat(rm),
                sp_mat(id),
                vec{5, 5},
                Col<char>{GRB_INTEGER, GRB_INTEGER, GRB_INTEGER, GRB_INTEGER},
                Col<char>{GRB_LESS_EQUAL, GRB_LESS_EQUAL}};

  Tree tree;
  tree.add_node(root);

  for (size_t s1 = 0; s1 != 21; ++s1)
  {
    for (size_t s2 = 0; s2 != 21; ++s2)
    {
      sub.d_rhs = {5.0 + 0.5*s1, 5.0 + 0.5*s2};
      tree.add_node(sub, 0);
    }
  }

  return tree;
}
