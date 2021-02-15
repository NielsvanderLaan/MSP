#include <iostream>

#include "main.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
  try
  {
    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_Threads, 1);

    int nstages = stoi(argv[1]);
    int n_outcomes = stoi(argv[2]);
    int depth = stoi(argv[3]);
    bool affine = stoi(argv[4]);

    cout << "number of stages: " << nstages <<
            "\noutcomes per stage: " << n_outcomes << '\n' <<
            (affine ? "Lagrangian" : "Scaled") << " cuts" <<
            "\ndepth: " << depth << '\n' << endl;

    Stagewise sw = ctrl_1D(env, nstages, n_outcomes);

    /*
    GRBModel sw_model = sw.lsde();
    sw_model.set(GRB_IntParam_OutputFlag, 1);
    sw_model.optimize();
    */

    cout << "SSDMIP\n";
    sw.decom(depth);
    sw.sddmip(affine);

  } catch (GRBException &e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << '\n';
    exit(e.getErrorCode());
  }
}
