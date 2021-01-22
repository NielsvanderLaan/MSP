#ifndef MSP_ENUMERATOR_H
#define MSP_ENUMERATOR_H

#include "gurobi_c++.h"
#include "../../nodedata/nodedata.h"
#include "../structs/structs.h"

using namespace std;

typedef vector<GRBVar> vvar;

class Enumerator
{
public:
  NodeData &d_data;
  GRBModel d_mp;
  GRBVar d_alpha;        // intercept
  vector<vvar> d_beta;   // d_beta[0] --> x_1, ..., x_a(n)
  vvar d_tau;            // d_tau[0] --> theta_1,..., theta_a(n)

  /*
   * objective of cgsp: c_n x_n + theta_n + beta.hat [x_a(n)] + tau.hat [theta_a(n)]
   */

  GRBModel d_sp;
  vector<vvar> d_x;      // x_1, ...., x_n
  vvar d_theta;          // theta_1, ..., theta_n

  vector<Solution> d_points;    // depth = n

  Enumerator(vector<NodeData> &nodes, vector<int> path, bool leaf, GRBEnv &env);
  Enumerator(Enumerator const &other);

  void add_cut(Cut &cut);
  void add_cut_to_sp(Cut &cut);
  void add_cut_to_mp(Cut &cut);

  Cut generate_cut(double rho, double tol = 1e-4);

  void solve_mp();
  void set_mp(Solution const &sol);
  void add_point(Solution point);    // c_n x_n + theta_n >= alpha - beta[x_a(n)] - tau[theta_a(n)]
  Cut candidate();

  void set_sub(Cut &cut);
  void solve_sp();
  Solution point();

  double crho();
  double alpha();
  double sub_val();
  double sub_bound();

  void set_rho(double rho);
};

#endif //MSP_ENUMERATOR_H
