#ifndef MSP_NODEDATA_H
#define MSP_NODEDATA_H

#include <armadillo>
#include "gurobi_c++.h"

using namespace std;

class NodeData
{
public:
  int d_stage;
  double d_prob;              // p_n (not the transition probability)
  double d_L;

  arma::vec d_costs;
  arma::vec d_lb;
  arma::vec d_ub;

  arma::sp_mat d_Amat;
  arma::sp_mat d_Bmat;
  arma::vec d_rhs;

  arma::Col<char> d_types;
  arma::Col<char> d_senses;

  GRBModel to_model(GRBEnv &env, bool lp = false);

  int nvars();
  int ncons();
};

#endif //MSP_NODEDATA_H
