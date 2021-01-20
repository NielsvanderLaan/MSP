#include <iostream>

#include "gurobi_c++.h"

#include "tree/tree.h"
#include "lsde/lsde.h"
#include "data/instances.h"

using namespace std;
using namespace arma;

int main()
{
  GRBEnv env;
  env.set(GRB_IntParam_OutputFlag, 0);

  Tree tree= ssv();

  tree.init_decom(env);

  Master &master = tree.d_masters[0];
  while (true)
  {
    master.optimize();

    Solution sol{vvec{master.xvals()}, vdouble{master.theta()}};

    Cut cut = tree.lp_cut(0, sol);

    if (not tree.add_cut(0, cut, sol))
      break;
  }

  cout << master.obj() << '\n';

  return 0;
}
