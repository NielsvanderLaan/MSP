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

  GRBModel d_sp;
  vector<vvar> d_x;      // x_1, ...., x_n
  vvar d_theta;          // theta_1, ..., theta_n

  vector<Solution> d_points;    // depth = n

  Enumerator(vector<NodeData> &nodes, vector<int> path, int mp_depth, bool leaf, GRBEnv &env);
  Enumerator(Enumerator const &other);

  void add_cut(Cut &cut);
  void add_cut_to_sp(Cut &cut);
  void add_cut_to_mp(Cut &cut);

  Cut opt_cut(double rho, double tol);
  Cut feas_cut(Solution const &sol, double tol);
  Cut fdecom(double tol);     // row generation / vertex enumeration
    // mp management
  void solve_mp();
  void set_mp(Solution const &sol);
  void add_point(Solution point);
  Cut candidate();
    // sp management
  void solve_sp();
  void set_sub(Cut &cut);
  Solution point();
    // getters
  double crho() const;
  double alpha() const;
  double sub_val() const;
  double sub_bound() const;
    // setters
  void set_rho(double rho);
  void set_bounds(double M = 1e3);
};

#endif //MSP_ENUMERATOR_H
