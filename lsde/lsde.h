#ifndef MSP_LSDE_H
#define MSP_LSDE_H

#include <queue>
#include <unordered_map>
#include <vector>
#include <armadillo>
#include "gurobi_c++.h"

#include "../tree/tree.h"

using namespace std;

class LSDE
{
public:
  GRBModel d_model;
  Tree &d_tree;

  unordered_map<int, vector<GRBVar>> d_vars;

  LSDE(GRBEnv &env, Tree &tree);
  void add(int id);

};

#endif //MSP_LSDE_H
