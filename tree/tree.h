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

  GRBModel lsde(GRBEnv &env);             // constructs lsde
  void decom(GRBEnv &env);                // initializes master problems etc.
  void solve();                           // solves using forward and backward passes
  void forward(bool lp);
  bool backward(int node = 0);

  bool solve_master(int node, bool lp = false);   // cutting plane approach to solve the master problem in node n
  bool add_cut(int id, Cut &cut);                 // adds cuts to master problem and propagates it through the tree.

    // auxiliary functions
  vvar add_to_lsde(int node, GRBModel &model, vvar const &parent_vars);


    //  optimality cuts
  Cut sddp_cut(int node);                             // standard Benders' cuts
  Cut scaled_cut(int node, double tol = 1e-4);        // scaled cuts (vertex enumeration)
  Cut cpt_scaled_cut(int node, double tol = 1e-4);    // scaled cuts (using Fenchel cutting planes)


  Cut fenchel_cut(int node, double tol = 1e-4);     // cutting plane for MIPs
  Cut fp_iteration(int node, vector<CutFamily> &gens, double tol = 1e-4);   // TODO

  void init_enums(int node);

  vector<int> path_to(int id);
  vector<int> nvars(int id);      // nvars of nodes on path to id
  bool is_leaf(int node);
};

#endif //MSP_TREE_H
