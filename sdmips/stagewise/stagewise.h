#ifndef MSP_STAGEWISE_H
#define MSP_STAGEWISE_H

#include "../nodedata/nodedata.h"
#include "../decomposition/master/master.h"
#include "../decomposition/enumerator/enumerator.h"
#include <random>

typedef vector<NodeData> stage_data;
typedef vector<Master> vmaster;
typedef vector<Enumerator> v_enum;
typedef vector<int> vpath;
typedef vector<Solution> vsol;

class Stagewise
{
public:
  mt19937 d_engine{random_device{}()};

  vector<stage_data> d_stages;
  vector<vmaster> d_masters;
  vector<v_enum> d_enumerators;
  vector<v_enum> d_fenchel;

  void add_node(NodeData const &data);

  GRBModel lsde(GRBEnv &env);
  void decom(GRBEnv &env);

  void sddmip();
  vector<vsol> forward(vector<vpath> &paths, bool lp);
  void backward(vector<vsol> const &sols);
  vector<vpath> sample(size_t nsamples = 30);
  vector<vpath> enumerate_paths(vector<vpath> paths = vector<vpath> (1));

  void solve(int stage, int node, bool lp, bool force);
  bool add_cp(Cut &cut, int stage, int node, double tol = 1e-4);
  void add_cut(Cut &cut, int stage);

  Cut sddp_cut(int stage, Solution const &sol);             // has to be valid for Q_t (stage-specific)
  Cut scaled_cut(int stage, Solution const &sol, double tol = 1e-4);
  Cut fenchel_cut(int stage, int node, double tol = 1e-4);                // has to be valid for X_n (node-specific)


  void init_enums(int stage, Solution const& sol);

  int outcomes(int stage) const;
  vector<int> nvars(int stage) const;
  vector<double> probs(stage_data const &stage) const;
};

#endif //MSP_STAGEWISE_H
