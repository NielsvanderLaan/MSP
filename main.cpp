#include <iostream>

#include "main.h"

using namespace std;
using namespace arma;

int main()
{
  GRBEnv env;
  env.set(GRB_IntParam_OutputFlag, 0);

  Tree tree = ssv();

  tree.init_decom(env);

  Master &master = tree.d_masters[0];
  /*
  while (true)
  {
    tree.forward(true);
    cout << master.obj() << '\n';
    Cut cut = tree.lp_cut(0);

    if (not tree.add_cut(0, cut))
      break;
  }
  cout << "sddp done\n";
  */
  while (true)
  {
    tree.forward(false);
    cout << master.obj() << '\n';

    if (not tree.backward())
      break;
  }

  return 0;
}
