#ifndef MSP_BENDERS_H
#define MSP_BENDERS_H

#include "../stagewise/stagewise.h"
#include "../decomposition/master/master.h"
#include "../decomposition/enumerator/enumerator.h"
#include "../decomposition/gomory/gomory.h"
#include "families.h"
#include <random>
#include <chrono>
#include <memory>

typedef vector<NodeData> stage_data;
typedef vector<Enumerator> v_enum;
typedef vector<int> vpath;
typedef vector<Solution> vsol;
typedef vector<Cut> outer_apx;



class Benders
{
public:
    GRBEnv &d_env;
    Stagewise &d_data;
    int d_depth;
    mt19937 d_engine;

    Benders(GRBEnv &env, Stagewise &data, int depth);
    virtual ~Benders() = default;

    void decom(Family type, size_t max_iter = 100, bool lp = false, size_t nsamples = 30, size_t eval = 1000);
    vector<vsol> forward(vector<vpath> const &paths, bool lp);
    void backward(Family type, vector<vsol> const &sols, vector<vpath> const &paths);

    vector<vpath> sample(size_t nsamples, size_t nstages = -1);
    vector<vpath> enumerate_paths(int start, int end, vector<vpath> const &paths = vector<vpath>(1));
    vector<vpath> enumerate_paths();

    Cut compute_cut(Family type, Solution const &sol, int stage, int node = 0);
    Cut sddp_cut(int stage, Solution const &sol);
    Cut scaled_cut(int stage, int node, Solution const &sol, bool affine, double tol = 1e-4);
    Cut lbda_cut(int stage, Solution const &sol, Alpha type, arma::vec alpha = {});

    v_enum &init_enums(int stage, int node, Solution const &sol);

    void add_cut(Cut &cut, Solution const &sol, bool shared, int stage, int node = -1);

    NodeData const &node_data(int stage, int outcome);
        // tree structure
    int master_idx(int stage, vpath const &path) const;
    int outcome(int stage, int node) const;
    vector<int> parents(int stage, int node) const;
    vector<int> children(int stage, int node) const;
    vector<int> tail(int stage, int node) const;

    double lb();
    double ub(size_t nsamples);

    virtual void add_cut(Cut &cut, int stage, int node) = 0;
    virtual void add_shared_cut(Cut &cut, int stage) = 0;

    virtual Master &get_master(int stage, int node) = 0;
    virtual v_enum &get_enums(int stage, int node) = 0;
    virtual Gomory &get_gomory(int stage, int out) = 0;

    void import_cuts(vector<outer_apx> stage_apx);
    virtual vector<outer_apx> export_cuts() = 0;
};

#endif //MSP_BENDERS_H
