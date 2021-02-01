#ifndef MSP_STAGEWISE_H
#define MSP_STAGEWISE_H

#include "../nodedata/nodedata.h"
#include "../decomposition/master/master.h"
#include "../decomposition/enumerator/enumerator.h"
#include <random>

typedef vector<NodeData> stage_data;
typedef vector<Master> vmaster;
typedef vector<int> path;
typedef vector<Solution> vsol;

class Stagewise
{
public:
  mt19937 d_engine{random_device{}()};

  vector<stage_data> d_stages;
  vector<vmaster> d_masters;

  void add_node(NodeData const &data);

  GRBModel lsde(GRBEnv &env);
  void decom(GRBEnv &env);

  void sddmip();
  vector<vsol> forward(vector<path> &paths);
  void backward(vector<vsol> const &sols);
  vector<path> sample(size_t nsamples = 30);
  vector<path> enumerate_paths(vector<path> paths = vector<path> (1));

  void solve(int stage, int node);
  void add_cut(Cut &cut, int stage);

  Cut sddp_cut(int stage, Solution const &sol);        // has to be valid for Q_t (stage-specific)
  Cut fenchel_cut(int stage, int node);                // has to be valid for X_n (node-specific)

  vector<int> nvars(int stage) const;
};

#endif //MSP_STAGEWISE_H
