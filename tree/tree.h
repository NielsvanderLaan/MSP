#ifndef MSP_TREE_H
#define MSP_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <assert.h>

#include "nodedata/nodedata.h"
#include "decomposition/master/master.h"

using namespace std;

class Tree
{

public:
    // use unordered maps?
  vector<NodeData> d_nodes;
  vector<Master> d_masters;
  vector<Enumerator> d_enumerators;

  unordered_map<int, int> d_ancestor;
  unordered_map<int, vector<int>> d_children;

  int add_node(NodeData &data, int ancestor = -1);    // ancestor = -1 --> root

  void init_decom(GRBEnv &env);

  vector<int> path_to(int id);
  vector<int> nvars(int id);      // nvars of nodes on path to id

  bool add_cut(int id, Cut &cut, Solution &sol);

  // cuts
  Cut lp_cut(int node, Solution &sol);


};

#endif //MSP_TREE_H
