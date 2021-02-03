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
  NodeData d_data;
  GRBModel *d_mp;
  GRBVar d_alpha;        // intercept
  vector<vvar> d_beta;   // d_beta[0] --> x_1, ..., x_a(n)
  vvar d_tau;            // d_tau[0] --> theta_1,..., theta_a(n)

  GRBModel *d_sp;
  vector<vvar> d_x;      // x_1, ...., x_n
  vvar d_theta;          // theta_1, ..., theta_n

  vector<Solution> d_points;    // depth = n
  vector<bool> d_directions;
  Solution d_prime;

  Enumerator(vector<NodeData> const &nodes, vector<int> path, size_t mp_depth, bool leaf, GRBEnv &env);
  Enumerator(Enumerator const &other);
  Enumerator(Enumerator &&other);
  ~Enumerator();

  void add_cut(Cut const &cut);
  void add_cut_to_sp(Cut const &cut);
  void add_cut_to_mp(Cut const &cut);

  Cut opt_cut(double rho, double tol);
  Cut feas_cut(Solution const &sol, double tol);
  Cut fdecom(double tol, bool reset = false);     // row generation / vertex enumeration
    // mp management
  void solve_mp();
  void clear();
  void set_mp(Solution const &sol);
  void add_point(Solution point, bool direction = false);
  Cut candidate();
    // sp management
  void solve_sp();
  void set_sub(Cut &cut);
  Solution point();
  Solution direction();
    // getters
  int sp_status() const;
  double crho() const;
  double alpha() const;
  double sub_val() const;
  double sub_bound() const;
  int mp_status() const;
  double mp_violation() const;
    // setters
  void set_rho(double rho);
  void set_bounds(double M = 1e3);
  void disable_tau();
  void prime(Solution const &point);
  void set_mp(bool tight);
};

#endif //MSP_ENUMERATOR_H
