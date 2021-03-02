#ifndef MSP_GOMORY_H
#define MSP_GOMORY_H

#include <memory>
#include "gurobi_c++.h"
#include "assert.h"

#include "../../nodedata/nodedata.h"
#include "../structs/structs.h"


using namespace std;

class Gomory
{
public:
  unique_ptr<GRBModel> d_model;
  NodeData const &d_data;
  GRBVar d_theta;
  vector<double> d_intercepts;

  Gomory(GRBEnv &env, NodeData const &data, bool leaf);
  Gomory(Gomory const &other);
  Gomory(Gomory &&other);

  void solve();
  double obj() const;   // lambda^T(omega - alpha) + psi(omega - alpha)
  void update(arma::vec rhs, vector<int> const &vbasis, vector<int> const &cbasis);
  void add_cut(Cut const &cut);
};

#endif //MSP_GOMORY_H
