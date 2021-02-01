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
  /*
  Tree tree = control_1D();
  GRBModel model = tree.lsde(env);
  model.set(GRB_IntParam_OutputFlag, 1);
  model.optimize();
  cout << "STOCHASTIC NESTED DECOMPOSITION\n";
  tree.decom(env);
  tree.SND();
  */

  Stagewise sw = ctrl_1D();
  GRBModel sw_model = sw.lsde(env);
  sw_model.set(GRB_IntParam_OutputFlag, 1);
  sw_model.optimize();
  cout << "SSDMIP\n";

  sw.decom(env);
  sw.sddmip();
}
