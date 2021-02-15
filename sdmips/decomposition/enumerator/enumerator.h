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
  GRBVar d_obj;
  GRBVar d_alpha;        // intercept
  vector<vvar> d_beta;   // d_beta[0] --> x_1, ..., x_a(n)
  vvar d_tau;            // d_tau[0] --> theta_1,..., theta_a(n)

  GRBConstr d_objcon;

  GRBModel *d_sp;
  vector<vvar> d_x;      // x_1, ...., x_n
  vvar d_theta;          // theta_1, ..., theta_n

  vector<Solution> d_points;    // depth = n

  Enumerator() = default;
  Enumerator(vector<NodeData> const &nodes, vector<int> path, size_t mp_depth, bool leaf, GRBEnv &env);
  Enumerator(vector<NodeData> const &nodes, vector<int> path, vector<Cut> const &cuts,
             size_t mp_depth, bool leaf, GRBEnv &env);
  Enumerator(Enumerator const &other);
  Enumerator(Enumerator &&other);
  ~Enumerator();

  void add_cut(Cut const &cut);
  void add_cut_to_sp(Cut const &cut);
  void add_cut_to_mp(Cut const &cut);

  Cut opt_cut(double rho, bool affine, double tol);
  Cut feas_cut(Solution const &sol, bool affine, double tol);
  Cut fdecom(double tol, bool affine, bool reset = false);
    // mp management
  void optimize_mp();
  Cut solve_mp(bool affine, double M);
  void clear();
  void set_mp(Solution const &sol);
  void add_point(Solution point, bool prime = false);
  Cut candidate();
    // sp management
  void solve_sp();
  void set_sub(Cut &cut);
  Solution point();
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
  void set_bounds(bool affine, double M);
  void set_tau_bounds(bool affine, double M);
  void set_mp(bool tight);
};

class mp_exception : public exception
{
  const char * what() const noexcept override
  {
    return "cgmp is numerically unstable\n";
  }
};

class sp_exception : public exception
{
  const char * what() const noexcept override
  {
    return "cgsp is numerically unstable\n";
  }
};


#endif //MSP_ENUMERATOR_H
