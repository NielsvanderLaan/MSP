#include "enumerator.h"

Enumerator::Enumerator(vector<NodeData> const &nodes, vector<int> path, size_t mp_depth, bool leaf, GRBEnv &env)
:
d_data(nodes[path.back()])
{
  d_mp = new GRBModel(env);
  d_sp = new GRBModel(env);
  d_sp->set(GRB_DoubleParam_MIPGapAbs, 1e-4);
  d_sp->set(GRB_DoubleParam_MIPGap, 0.0);
  d_sp->set(GRB_DoubleParam_TimeLimit, 60);
  d_mp->set(GRB_IntParam_Method, 1);
  d_mp->set(GRB_IntAttr_ModelSense, -1);
  d_mp->set(GRB_DoubleParam_TimeLimit, 10);

  d_obj = d_mp->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
  d_alpha = d_mp->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

  d_objcon = d_mp->addConstr(d_obj == d_alpha, "obj");

  for (size_t stage = 0; stage < mp_depth; ++stage)
  {
    int nvars = nodes[path[stage]].nvars();
    vector<double> lb(nvars, -GRB_INFINITY);
    GRBVar *beta  = d_mp->addVars(lb.data(),
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nvars);
    d_beta.push_back(vvar(beta, beta + nvars));
    delete[] beta;
  }

  for (size_t stage = 0; stage < mp_depth; ++stage)
    d_tau.push_back(d_mp->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "tau_" + to_string(stage)));


  for (int node : path)
  {
    NodeData const &data = nodes[node];
    GRBVar *xnode = d_sp->addVars(data.d_lb.memptr(),
                                  data.d_ub.memptr(),
                                  data.d_costs.memptr(),
                                  data.d_types.memptr(),
                                 nullptr,
                                  data.nvars());
    d_x.push_back(vvar(xnode, xnode + data.nvars()));
    delete[] xnode;
  }

  for (auto it = path.begin(); it != path.end(); ++it)
    d_theta.push_back(d_sp->addVar(nodes[*it].d_L,
                                  GRB_INFINITY,
                                  1.0,
                                  GRB_CONTINUOUS));

  if (leaf)
  {
    d_theta.back().set(GRB_DoubleAttr_LB, 0);
    d_theta.back().set(GRB_DoubleAttr_UB, 0);
  }

  for (size_t stage = 0; stage != path.size(); ++stage)
  {
    NodeData const &data = nodes[path[stage]];
    GRBLinExpr lhs[data.ncons()];

    for (auto iter = data.d_Amat.begin(); iter != data.d_Amat.end(); ++iter)
      lhs[iter.col()] += *iter * d_x[stage][iter.row()];

    if (stage > 0)   // check if not root (ancestor variables have to exist)
    {
      for (auto iter = data.d_Bmat.begin(); iter != data.d_Bmat.end(); ++iter)
        lhs[iter.col()] += *iter * d_x[stage - 1][iter.row()];
    }

    delete[] d_sp->addConstrs(lhs,
                             data.d_senses.memptr(),
                             data.d_rhs.memptr(),
                             nullptr,
                             data.ncons());
  }
  d_mp->update();
  d_sp->update();
}

Enumerator::Enumerator(const Enumerator &other)
:
d_data(other.d_data),
d_points(other.d_points)
{
  d_mp = new GRBModel(*other.d_mp);
  d_sp = new GRBModel(*other.d_sp);

  GRBVar *mp_vars = d_mp->getVars();
  d_obj = mp_vars[0];
  d_alpha = mp_vars[1];
  int start = 2;
  for (auto it = other.d_beta.begin(); it != other.d_beta.end(); ++it)
  {
    d_beta.push_back(vvar(mp_vars + start,
                          mp_vars + start + it->size()));
    start += it->size();
  }

  for (size_t var = 0; var != other.d_tau.size(); ++var)
    d_tau.push_back(mp_vars[start + var]);

  d_objcon = d_mp->getConstrByName("obj");

  GRBVar *sub_vars = d_sp->getVars();

  start = 0;
  for (auto it = other.d_x.begin(); it != other.d_x.end(); ++it)
  {
    d_x.push_back(vvar(sub_vars + start,
                       sub_vars + start + it->size()));
    start += it->size();
  }

  for (size_t var = 0; var != other.d_theta.size(); ++var)
    d_theta.push_back(sub_vars[start + var]);

  delete[] mp_vars;
  delete[] sub_vars;

  d_mp->update();
  d_sp->update();
}

Enumerator::Enumerator(Enumerator &&other)
:
d_data(other.d_data),
d_mp(other.d_mp),
d_obj(other.d_obj),
d_alpha(other.d_alpha),
d_beta(other.d_beta),
d_tau(other.d_tau),
d_objcon(other.d_objcon),
d_sp(other.d_sp),
d_x(other.d_x),
d_theta(other.d_theta),
d_points(other.d_points)
{
  other.d_sp = nullptr;
  other.d_mp = nullptr;
}

Enumerator::~Enumerator()
{
  if (d_mp)
    delete d_mp;

  if (d_sp)
    delete d_sp;
}








