#include "enumerator.h"

Enumerator::Enumerator(vector<NodeData> &nodes, vector<int> path, bool leaf, GRBEnv &env)
:
d_mp(env),
d_sp(env)
{
  if (path.size() == 1)   // root node
    return;

  d_alpha = d_mp.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

  for (auto it = path.begin() + 1; it != path.end(); ++it)
  {
    int nvars = nodes[*it].nvars();
    vector<double> lb(nvars, -GRB_INFINITY);
    GRBVar *beta  = d_mp.addVars(lb.data(),
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nvars);
    d_beta.push_back(vvar(beta, beta + nvars));
    delete[] beta;
  }

  for (size_t idx = 1; idx != path.size(); ++idx)
    d_tau.push_back(d_mp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));

  for (int node : path)
  {
    int nvars = nodes[node].nvars();
    GRBVar *xnode = d_sp.addVars(nodes[node].d_lb.memptr(),
                                 nodes[node].d_lb.memptr(),
                                 nodes[node].d_costs.memptr(),     // only relevant for x_n
                                 nodes[node].d_types.memptr(),
                                 nullptr,
                                 nvars);
    d_x.push_back(vvar(xnode, xnode + nvars));
    delete[] xnode;
  }

  for (auto it = path.begin() + leaf; it != path.end(); ++it)      // first is skipped if leaf = true{
    d_theta.push_back(d_sp.addVar(nodes[*it].d_L,
                                  GRB_INFINITY,
                                  1.0,                         // only relevant for theta_n
                                  GRB_CONTINUOUS));



  for (size_t lvl = 0; lvl != path.size(); ++lvl)
  {
    NodeData &data = nodes[path[lvl]];
    GRBLinExpr lhs[data.ncons()];

    for (auto iter = data.d_Amat.begin(); iter != data.d_Amat.end(); ++iter)
      lhs[iter.col()] += *iter * d_x[lvl][iter.row()];

    if (lvl + 1 != path.size())   // check if not root (ancestor variables have to exist)
    {
      for (auto iter = data.d_Bmat.begin(); iter != data.d_Bmat.end(); ++iter)
        lhs[iter.col()] += *iter * d_x[lvl + 1][iter.row()];
    }

    delete[] d_sp.addConstrs(lhs,
                             data.d_senses.memptr(),
                             data.d_rhs.memptr(),
                             nullptr,
                             data.ncons());
  }

  d_mp.update();
  d_sp.update();
}

Enumerator::Enumerator(const Enumerator &other)
:
d_mp(other.d_mp),
d_sp(other.d_sp),
d_points(other.d_points)
{
  if (d_mp.get(GRB_IntAttr_NumVars) == 0)   // root node
    return;

  GRBVar *mp_vars = d_mp.getVars();
  d_alpha = mp_vars[0];
  int start = 1;
  for (auto it = other.d_beta.begin(); it != other.d_beta.end(); ++it)
  {
    d_beta.push_back(vvar(mp_vars + start,
                          mp_vars + start + it->size()));
    start += it->size();
  }

  for (size_t var = 0; var != other.d_tau.size(); ++var)
    d_tau.push_back(mp_vars[start + var]);
  GRBVar *sub_vars = d_sp.getVars();

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

  d_mp.update();
  d_sp.update();
}












