#include <iostream>

#include "main.h"

using namespace std;


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
    bool sparse = (argc > 5 && string(argv[5]) == "SPARSE");

    cout << "number of stages: " << nstages <<
            "\noutcomes per stage: " << n_outcomes << '\n' <<
            (affine ? "Lagrangian" : "Scaled") << " cuts" <<
            "\ndepth: " << depth << '\n' <<
            (sparse ? "sparse" : "dense") << '\n' << endl;

    Stagewise sw = ctrl_1D(nstages, n_outcomes);
    //Stagewise sw = sclsp(nstages, n_outcomes);

    /*
    GRBModel sw_model = sw.lsde(env);
    sw_model.set(GRB_IntParam_OutputFlag, 1);
    sw_model.optimize();
    */

    unique_ptr<Benders> benders;
    if (sparse)
      benders = make_unique<spBenders>(env, sw, depth);
    else
      benders = make_unique<dBenders>(env, sw, depth, 10);

    cout << "SDDP" << endl;
    benders->decom(SDDP, 5, false);
    cout << "SSDMIP" << endl;
    benders->decom(affine ? LR : SC, 250, false);

  } catch (GRBException &e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << '\n';
    exit(e.getErrorCode());
  }
}
