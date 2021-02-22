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

    print(argc, argv);
    Args args = parse(argc, argv);
    Stagewise sw = get_problem(args);

    run(env,sw, args.ws_types, args.types, args.iter_limits, args.depth, args.sparse, args.samples);

    /*
    GRBModel sw_model = sw.lsde(env);
    sw_model.set(GRB_IntParam_OutputFlag, 1);
    sw_model.optimize();
    */

  } catch (GRBException &e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << '\n';
    exit(e.getErrorCode());
  }
}
