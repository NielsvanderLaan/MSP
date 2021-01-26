#ifndef MSP_TREE_H
#define MSP_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <assert.h>

#include "nodedata/nodedata.h"

#include "decomposition/cutfamily/cutfamily.h"
#include "decomposition/master/master.h"
#include "decomposition/enumerator/enumerator.h"

using namespace std;

class Tree
{

public:
  vector<NodeData> d_nodes;
  vector<Master> d_masters;
  vector<Enumerator> d_enumerators;
  vector<Enumerator> d_fenchel;

  unordered_map<int, int> d_ancestor;
  unordered_map<int, vector<int>> d_children;

      // for constructing the problem
  int add_node(NodeData &data, int ancestor = -1);    // ancestor = -1 --> root

    // constructs the LSDE
  GRBModel lsde(GRBEnv &env);

  void init_decom(GRBEnv &env);
  void forward(bool lp);
  bool backward(int node = 0);

  void solve_master(int node, bool lp);


    // auxiliary functions
  vvar add_to_lsde(int node, GRBModel &model, vvar const &parent_vars);

  vector<int> path_to(int id);
  vector<int> nvars(int id);      // nvars of nodes on path to id

  bool add_cut(int id, Cut &cut);

  // cuts
  Cut lp_cut(int node);   // TODO: implement fixed point iteration
  Cut scaled_cut(int node, double tol = 1e-4);
  Cut fenchel_cut(int node, double tol = 1e-4);

  Cut fp_iteration(int node, vector<CutFamily> &gens, double tol = 1e-4);

  void init_enums(int node);

  bool is_leaf(int node);
};

#endif //MSP_TREE_H
