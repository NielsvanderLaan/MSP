#include "gomory.h"

Gomory::Gomory(GRBEnv &env, NodeData const &data, bool leaf)
:
d_model(make_unique<GRBModel>(data.to_model(env))),
d_data(data)
{
  d_theta = d_model->addVar(leaf ? 0.0 : d_data.d_L,
                            leaf? 0.0 : GRB_INFINITY,
                            1.0,
                            GRB_CONTINUOUS,
                            "theta");
  d_model->set(GRB_DoubleParam_TimeLimit, 1.0);
  d_model->update();
}

Gomory::Gomory(Gomory const &other)
:
d_model(make_unique<GRBModel>(*other.d_model)),
d_data(other.d_data),
d_intercepts(other.d_intercepts)
{
  d_theta = d_model->getVarByName("theta");
}

Gomory::Gomory(Gomory &&other)
:
d_model(move(other.d_model)),
d_data(other.d_data),
d_theta(other.d_theta),
d_intercepts(move(other.d_intercepts))
{}

void Gomory::solve()
{
  d_model->optimize();
  int status = d_model->get(GRB_IntAttr_Status);
  assert(status == GRB_OPTIMAL);
}

void Gomory::add_cut(Cut const &cut)
{
  d_intercepts.push_back(cut.d_alpha);

  GRBVar *xvars = d_model->getVars();
  int stage = d_data.d_stage - 1;
  GRBLinExpr lhs = (1 + cut.d_tau[stage]) * d_theta;
  lhs.addTerms(cut.d_beta[stage].memptr(), xvars, d_data.nvars());

  d_model->addConstr(lhs >= cut.d_alpha);
  delete[] xvars;
  d_model->update();
}

void Gomory::update(arma::vec rhs, vector<int> const &basis)
{
  arma::vec lb = d_data.d_lb;
  int ncons = rhs.size();
  arma::vec cut_rhs(d_intercepts);

  d_theta.set(GRB_DoubleAttr_LB, d_data.d_L);

  for (int basic : basis)
  {
    if (basic < lb.size())
    {
      lb[basic] = -GRB_INFINITY;
      continue;
    }

    basic -= lb.size();
    if (basic == 0)     // basic corresponds to theta
    {
      d_theta.set(GRB_DoubleAttr_LB, -GRB_INFINITY);
      continue;
    }

    --basic;
    if (basic < ncons)     // basic corresponds to a slack variable of a constraint
    {
      switch (d_data.d_senses[basic])
      {
        case GRB_EQUAL:
          break;
        case GRB_LESS_EQUAL:
          rhs[basic] = GRB_INFINITY;
          break;
        case GRB_GREATER_EQUAL:
          rhs[basic] = -GRB_INFINITY;
      }
      continue;
    }

    basic -= ncons;
    cut_rhs[basic] = -GRB_INFINITY;
  }

  GRBVar* vars = d_model->getVars();
  d_model->set(GRB_DoubleAttr_LB, vars, lb.memptr(), lb.n_elem);
  delete[] vars;

  GRBConstr* cons = d_model->getConstrs();
  d_model->set(GRB_DoubleAttr_RHS, cons, rhs.memptr(), ncons);
  d_model->set(GRB_DoubleAttr_RHS, cons + ncons, cut_rhs.memptr(), cut_rhs.size());
  delete[] cons;
}

void Gomory::update(arma::vec rhs, vector<int> const &vbasis, vector<int> const &cbasis)
{
  int ncons = d_data.ncons();
  for (size_t con = 0; con != ncons; ++con)
  {
    if (cbasis[con])
      continue;

    switch (d_data.d_senses[con])
    {
      case GRB_EQUAL:
        continue;
      case GRB_LESS_EQUAL:
        rhs[con] = GRB_INFINITY;
        break;
      case GRB_GREATER_EQUAL:
        rhs[con] = -GRB_INFINITY;
    }
  }

  arma::vec cut_rhs(d_intercepts);
  for (size_t cut = 0; cut != cut_rhs.size(); ++cut)
  {
    if (not cbasis[ncons + cut])
      cut_rhs[cut] = -GRB_INFINITY;
  }

  GRBConstr* cons = d_model->getConstrs();
  d_model->set(GRB_DoubleAttr_RHS, cons, rhs.memptr(), ncons);
  d_model->set(GRB_DoubleAttr_RHS, cons + ncons, cut_rhs.memptr(), cut_rhs.size());
  delete[] cons;


  arma::vec lb = d_data.d_lb;
  arma::vec ub = d_data.d_ub;

  for (size_t var = 0; var != d_data.nvars(); ++var)
  {
    if (vbasis[var] != -1)
      lb[var] = -GRB_INFINITY;
    if (vbasis[var] != -2)
      ub[var] = GRB_INFINITY;
  }

  GRBVar* vars = d_model->getVars();
  d_model->set(GRB_DoubleAttr_LB, vars, lb.memptr(), lb.n_elem);
  d_model->set(GRB_DoubleAttr_UB, vars, ub.memptr(), lb.n_elem);
  delete[] vars;

  d_theta.set(GRB_DoubleAttr_LB, vbasis.back() == -1 ? d_data.d_L : -GRB_INFINITY);

  d_model->update();
}

double Gomory::obj() const
{
  return d_model->get(GRB_DoubleAttr_ObjBound);
}
