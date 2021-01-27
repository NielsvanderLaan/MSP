#include "instances.h"

Tree control_1D()
{
  vector<int> scenarios {1, 10, 10}; // per stage
  int stages = scenarios.size();
  vector<vdouble> nodes(stages);       // stores the nodes for each stage

  sp_mat Amat(mat{{1,0}, {-1, 0}, {-1, 1}, {1, 1}});
  vec lb {0.0,           0.0,         0.0, 0.0};
  vec ub {GRB_INFINITY, GRB_INFINITY, 1.0, 1.0};
  Col<char> types {GRB_CONTINUOUS, GRB_CONTINUOUS, GRB_CONTINUOUS, GRB_CONTINUOUS};
  Col<char> senses {GRB_EQUAL, GRB_EQUAL};

  NodeData root {1,
                 1,
                 0,
                 vec{1.0, 1.0, 0.0, 0.0},
                 lb,
                 ub,
                 Amat,
                 sp_mat(),
                 vec{2.0, 1.0},
                 types,
                 senses };

  Tree tree;
  nodes[0].push_back(tree.add_node(root));

  sp_mat Bmat = {umat{{0, 1}, {0, 0}}, vec{-1.0, 1.0}, 4, 2};

  double beta = 0.9;
  double prob = 1.0;
  for (int stage = 1; stage != stages ; ++stage)
  {
    prob /= scenarios[stage];

    NodeData sub {stage + 1,
                  prob,
                  0,
                  vec{beta, beta, 0.0, 0.0},
                  vec{0.0,           0.0,         0.0, 0.0},
                  vec{GRB_INFINITY, GRB_INFINITY, 1.0, 1.0},
                  Amat,
                  Bmat,
                  vec{0.0, 1.0},
                  types,
                  senses};

    for (int parent : nodes[stage - 1])
    {
      for (int s = 0; s != scenarios[stage]; ++s)
      {
        sub.d_rhs[0] = -4.5 + s;
        nodes[stage].push_back(tree.add_node(sub, parent));
      }
    }
    beta = beta * beta;
  }

  return tree;
}
