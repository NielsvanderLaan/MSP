#ifndef MSP_NODEDATA_H
#define MSP_NODEDATA_H

#include <vector>
#include <armadillo>
#include "gurobi_c++.h"
#include "assert.h"

using namespace std;

typedef vector<GRBVar> vvar;

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

  arma::uvec d_fixed_constrs;   // deterministic constraints
  arma::uvec d_box_constrs;     // deterministic non-linking constraints (only x_t, no x_{t - 1})
  size_t d_rn_constrs;

  GRBModel to_model(GRBEnv &env, bool lp = false) const;
  void to_box(bool wide);
  void bnds_to_cons();
  void clear();

  vvar add_to_lsde(GRBModel &lsde, vvar const& parent_vars, double corr = 1.0) const;

  int nvars() const;
  int ncons() const;
};

#endif //MSP_NODEDATA_H
