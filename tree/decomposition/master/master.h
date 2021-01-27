#ifndef MSP_MASTER_H
#define MSP_MASTER_H

#include "gurobi_c++.h"
#include "assert.h"

#include "../../nodedata/nodedata.h"
#include "../structs/structs.h"

typedef vector<GRBVar> vvar;

using namespace std;

class Master
{
public:
  NodeData &d_data;

  GRBModel d_mip;
  GRBVar d_mip_theta;
  vvar d_mip_xvars;

  GRBModel d_lp;            // nodal relaxation
  GRBVar d_lp_theta;        // costs-to-go
  vvar d_lp_xvars;          // decision variables

  vector<Cut> d_cuts;
  Solution d_state;         // state variables [x_a(n)], [theta_a(n)]

  Master(NodeData &data, bool leaf, GRBEnv &env);
  Master(const Master &other);

  Cut opt_cut();

  bool add_cut(Cut const &cut, double tol = 1e-4);
  void update(Solution const &sol);
  void solve_lp();
  void solve_mip();

  void set_rho(double rho);

    // getters
  arma::vec lp_xvals();
  double lp_theta();
  Solution lp_forward();
  bool integer();

  arma::vec mip_xvals();
  double mip_theta();
  Solution mip_forward();

  arma::vec multipliers();
  double obj();
};

#endif //MSP_MASTER_H
