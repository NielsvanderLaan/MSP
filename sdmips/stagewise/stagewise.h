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

typedef pair<Master, v_enum*> node;
typedef vector<node> vnode;


class Stagewise
{
public:
  mt19937 d_engine{891}; // random_device{}()

  vector<stage_data> d_stages;
  vector<vnode> d_nodes;
  int d_depth;

  void add_node(NodeData const &data);

  ~Stagewise();

  GRBModel lsde(GRBEnv &env);
  void decom(GRBEnv &env, int depth);

  void sddmip(bool affine);
  void sddp();                // TODO
  vector<vsol> forward(vector<vpath> const &paths, bool affine, bool lp);
  void backward(vector<vsol> const &sols, vector<vpath> const &paths, bool affine);
  void shared_backward(vector<vsol> const &sols, bool affine);

  vector<vpath> sample(size_t nsamples = 30);
  vector<vpath> enumerate_paths(int start, int end, vector<vpath> paths = vector<vpath>(1));
  vector<vpath> enumerate_paths();

  void solve_lp(int stage, int node);
  void solve_mip(int stage, int node);
  Solution solution(int stage, int node);
  void add_cut(Cut &cut, int stage, vector<int> const &path);
  void add_shared_cut(Cut &cut, int stage);

  Cut sddp_cut(int stage, Solution const &sol);         // sddp cuts are shared cuts by default
  Cut scaled_cut(int stage, int node, Solution const &sol, bool affine, double tol = 1e-4);
  Cut shared_scaled_cut(int stage, Solution const &sol, bool affine, double tol = 1e-4);

  void init_enums(int stage, int node, Solution const &sol);

  int master_idx(int stage, vpath const &path);
  vector<int> parents(int stage, vector<int> const &path);
  vector<int> children(int stage, int node);
  Master& get_master(int stage, int node);
  v_enum &get_enums(int stage, int node);
  int outcomes(int stage) const;
  vector<int> nvars(int stage) const;
  vector<double> probs(stage_data const &stage) const;
};

#endif //MSP_STAGEWISE_H
