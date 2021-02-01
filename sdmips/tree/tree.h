#ifndef MSP_TREE_H
#define MSP_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <assert.h>

#include "../stagewise/stagewise.h"
#include "../nodedata/nodedata.h"

#include "../decomposition/master/master.h"
#include "../decomposition/enumerator/enumerator.h"

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

  Tree() = default;
  Tree (Stagewise const &sw);
      // for constructing the problem
  int add_node(NodeData const &data, int ancestor = -1);    // ancestor = -1 --> root

  GRBModel lsde(GRBEnv &env);             // constructs lsde
  void decom(GRBEnv &env);                // initializes master problems etc.
  void SND();                           // solves using forward and backward passes
  void forward(bool lp);
  bool backward(int node = 0);

  void solve_master(int node, bool force, bool lp);   // cutting plane approach to SND the master problem in node n
  bool add_cut(int id, Cut &cut, double tol = 1e-4);  // adds cuts to master problem and propagates it through the sdmips.

    //  optimality cuts
  Cut sddp_cut(int node);                             // standard Benders' cuts
  Cut scaled_cut(int node, double tol = 1e-4);        // scaled cuts (vertex enumeration)
  Cut cpt_scaled_cut(int node, double tol = 1e-4);    // scaled cuts (using Fenchel cutting planes)
    // cutting planes for MIPs
  Cut fenchel_cut(int node, double tol = 1e-4);

  void init_enums(int node);

  vector<int> path_to(int id);
  vector<int> nvars(int id);      // nvars of nodes on path to id
  bool is_leaf(int node);
};

#endif //MSP_TREE_H
