#include "instances.h"

Stagewise ctrl_1D(size_t nstages, size_t n_outcomes)
{
  int seed = 1234; //random_device{}()
  mt19937 engine(seed);
  uniform_real_distribution<double> uni(0.5, 1.0);

  vector<int> scenarios {1};
  for (size_t stage = 0; stage != nstages - 1; ++stage)
    scenarios.push_back(n_outcomes);

  int stages = scenarios.size();

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


  Stagewise sw;
  sw.add_node(root);

  sp_mat Bmat = {umat{{0, 1}, {0, 0}}, vec{-1.0, 1.0}, 4, 2};

  double beta = 0.9;
  double costs = beta;
  for (int stage = 1; stage != stages ; ++stage)
  {
    double prob = 1.0 / scenarios[stage];

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
                  senses,
                  uvec{1},
                  uvec{1}};

    //double step = 10.0 / (scenarios[stage] - 1);
    for (int s = 0; s != scenarios[stage]; ++s)
    {
      //sub.d_rhs[0] = -5.0 + s * step;
      if (s % 2 == 0)
        sub.d_rhs[0] = uni(engine);
      else
        sub.d_rhs[0] = -uni(engine);
      sw.add_node(sub);
    }

    costs *= beta;
  }

  return sw;
}
