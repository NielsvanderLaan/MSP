#ifndef MSP_STAGEWISE_H
#define MSP_STAGEWISE_H

#include "../nodedata/nodedata.h"
#include "../decomposition/master/master.h"
#include "../decomposition/enumerator/enumerator.h"
#include <tuple>
#include <random>

typedef vector<NodeData> stage_data;
typedef vector<Master> vmaster;
typedef vector<Enumerator> v_enum;
typedef vector<int> vpath;
typedef vector<Solution> vsol;

typedef tuple<Master, Enumerator, v_enum*> node;
typedef vector<node> vnode;


class Stagewise
{
public:
  mt19937 d_engine{891}; // random_device{}()

  vector<stage_data> d_stages;
  vector<vnode> d_nodes;
  int d_depth;

  ~Stagewise();

  void add_node(NodeData const &data);

  GRBModel lsde(GRBEnv &env);
  void decom(GRBEnv &env, int depth);
  void decom(GRBEnv &env);            // stagewise outer approximations

  void sddmip(bool affine);
  void sddp();
  vector<vsol> forward(vector<vpath> const &paths, bool affine, bool lp);
  void backward(vector<vsol> const &sols, vector<vpath> const &paths, bool affine);
  void shared_backward(vector<vsol> const &sols, bool affine);

  vector<vpath> sample(size_t nsamples = 30);
  vector<vpath> enumerate_paths(int start, int end, vector<vpath> paths = vector<vpath>(1));
  vector<vpath> enumerate_paths();

  void solve(int stage, int node, bool affine, bool lp, bool force);
  Solution solution(int stage, int node);
  bool add_cp(Cut &cut, int stage, int node, double tol = 1e-4);
  void add_cut(Cut &cut, int stage, vector<int> const &path);
  void add_shared_cut(Cut &cut, int stage);

  Cut sddp_cut(int stage, Solution const &sol);         // sddp cuts are shared cuts by default
  Cut scaled_cut(int stage, int node, Solution const &sol, bool affine, double tol = 1e-4);
  Cut shared_scaled_cut(int stage, Solution const &sol, bool affine, double tol = 1e-4);
  Cut fenchel_cut(int stage, int node, bool affine, double tol = 1e-4);     // has to be valid for X_n (node-specific)

  void init_enums(int stage, int node, Solution const &sol);

  int master_idx(int stage, vpath const &path);
  vector<int> parents(int stage, vector<int> const &path);
  vector<int> children(int stage, int node);
  Master& get_master(int stage, int node);
  Enumerator &get_fenchel(int stage, int node);
  v_enum &get_enums(int stage, int node);
  int outcomes(int stage) const;
  vector<int> nvars(int stage) const;
  vector<double> probs(stage_data const &stage) const;
};

#endif //MSP_STAGEWISE_H
