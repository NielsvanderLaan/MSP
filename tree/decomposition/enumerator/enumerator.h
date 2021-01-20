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
  GRBModel d_mp;
  GRBVar d_alpha;        // intercept
  vector<vvar> d_beta;   // d_beta[0] --> x_a(n), ..., ... --> x_1
  vvar d_tau;            // d_tau[0] --> theta_a(n), ..., ... --> theta_1

  /*
   * objective of cgsp: c_n x_n + theta_n + beta.hat [x_a(n)] + tau.hat [theta_a(n)]
   */

  GRBModel d_sp;
  vector<vvar> d_x;      // x_n, ...., x_1
  vvar d_theta;          // theta_n, ..., theta_1

  /*
   * depth of points is n
   */

  vector<Solution> d_points;

  Enumerator(vector<NodeData> &nodes, vector<int> path, bool leaf, GRBEnv &env);
  Enumerator(Enumerator const &other);



  // the feasible region of the CGSP depends on the outer approximations of the ancestors
  // its objective is in terms of v.hat, which depends on the outer approximation in the current node
  // thus if a cut is added to one of the master problems, it should be propagated to its children (recursively)



  void add_cut(Cut &cut);
  void add_cut_to_sp(Cut &cut);
  void add_cut_to_mp(Cut &cut);

  void generate_cut();

  void solve_mp();
  void add_point(Solution &point);    // c_n x_n + theta_n >= alpha - beta[x_a(n)] - tau[theta_a(n)]
  Cut candidate();

  void set_sub(Cut &cut);
  void solve_sub();
  Solution point();

  double crho();
};

#endif //MSP_ENUMERATOR_H
