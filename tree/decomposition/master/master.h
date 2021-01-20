#ifndef MSP_MASTER_H
#define MSP_MASTER_H

#include "gurobi_c++.h"
#include "assert.h"

#include "../../nodedata/nodedata.h"
#include "../structs/structs.h"
#include "../enumerator/enumerator.h"

using namespace std;

class Master
{
public:
  NodeData &d_data;

  GRBModel d_model;         // nodal relaxation
  GRBVar d_theta;           // costs-to-go
  vector<GRBVar> d_xvars;   // decision variables

  vector<Cut> d_cuts;

  Master(NodeData &data, bool leaf, GRBEnv &env);
  Master(const Master &other);

  bool add_cut(Cut &cut, Solution &sol, double tol = 1e-4);
  void update(Solution &sol);
  void optimize();
  Cut compute_cut(Solution &sol);

    // getters
  arma::vec xvals();
  double theta();
  arma::vec multipliers();
  double obj();
};

#endif //MSP_MASTER_H
