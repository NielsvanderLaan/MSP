#include <iostream>

#include "main.h"

using namespace std;
using namespace arma;

int main()
{
  GRBEnv env;
  env.set(GRB_IntParam_OutputFlag, 0);
  env.set(GRB_IntParam_Threads, 1);

  //Tree tree = ssv();
  Tree tree = control_1D();
    /*
  GRBModel model = tree.lsde(env);
  model.set(GRB_IntParam_OutputFlag, 1);
  model.optimize();
    */
  cout << "running sddp\n";
  tree.decom(env);
  tree.solve();
}
