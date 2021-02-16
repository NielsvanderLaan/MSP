#include "instances.h"

Tree control_1D()
{
  int seed = 1234; //random_device{}()
  mt19937 engine(seed);
  uniform_real_distribution<double> uni(0.5, 1.0);

  vector<int> scenarios {1, 4, 4, 4, 4};   // per stage
  int stages = scenarios.size();
  vector<vdouble> nodes(stages);       // stores the nodes for each stage

  sp_mat Amat(mat{{1,0}, {-1, 0}, {-1, 1}, {1, 1}});
  double M = 10.0;
  vec lb {0, 0, 0, 0};
  vec ub {M, M, 1, 1};
  Col<char> types {GRB_CONTINUOUS, GRB_CONTINUOUS, GRB_INTEGER, GRB_INTEGER};
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
  double costs = beta;
  double prob = 1.0;
  for (int stage = 1; stage != stages ; ++stage)
  {
    prob /= scenarios[stage];

    NodeData sub {stage + 1,
                  prob,
                  0,
                  vec{costs, costs, 0.0, 0.0},
                  lb,
                  ub,
                  Amat,
                  Bmat,
                  vec{0.0, 1.0},
                  types,
                  senses};

    for (int s = 0; s != scenarios[stage]; ++s)
    {
      if (s % 2 == 0)
        sub.d_rhs[0] = uni(engine);
      else
        sub.d_rhs[0] = -uni(engine);

      for (int parent : nodes[stage - 1])
        nodes[stage].push_back(tree.add_node(sub, parent));

    }
    costs *= beta;
  }

  return tree;
}
