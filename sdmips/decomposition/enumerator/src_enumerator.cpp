#include "enumerator.h"

void Enumerator::solve_mp()
{
  d_mp.optimize();
  assert(d_mp.get(GRB_IntAttr_Status) == 2);
}

void Enumerator::solve_sp()
{
  d_sp.optimize();
}

int Enumerator::sp_status() const
{
  return d_sp.get(GRB_IntAttr_Status);
}

double Enumerator::sub_val() const
{
  return d_sp.get(GRB_DoubleAttr_ObjVal);
}

double Enumerator::sub_bound() const
{
  return d_sp.get(GRB_DoubleAttr_ObjBound);
}

double Enumerator::alpha() const
{
  return d_alpha.get(GRB_DoubleAttr_X);
}

double Enumerator::crho() const
{
  return d_mp.get(GRB_DoubleAttr_ObjVal) + alpha() - sub_bound();
}

void Enumerator::set_rho(double rho)
{
  d_tau.back().set(GRB_DoubleAttr_Obj, rho);
}

void Enumerator::disable_tau()
{
  vector<double> zeros(d_tau.size(), 0.0);
  d_mp.set(GRB_DoubleAttr_UB, d_tau.data(), zeros.data(), zeros.size());
}

void Enumerator::set_mp(Solution const &sol)
{
  int depth = sol.depth();
  for (size_t stage = 0; stage != depth; ++stage)
    d_mp.set(GRB_DoubleAttr_Obj,
             d_beta[stage].data(),
             sol.d_x[stage].memptr(),
             d_beta[stage].size());

  d_mp.set(GRB_DoubleAttr_Obj,
           d_tau.data(),
           sol.d_theta.data(),
           depth);

  d_mp.update();
}

Cut Enumerator::candidate()
{
  vvec beta_vals;
  beta_vals.reserve(d_beta.size());

  for (auto &beta : d_beta)
  {
    double *vals = d_mp.get(GRB_DoubleAttr_X, beta.data(), beta.size());
    beta_vals.push_back(arma::vec(vals, beta.size()));
    delete[] vals;
  }

  double *vals = d_mp.get(GRB_DoubleAttr_X, d_tau.data(), d_tau.size());
  vdouble tau_vals(vals, vals + d_tau.size());
  delete[] vals;

  return Cut {alpha(),
              beta_vals,
              tau_vals};
}

void Enumerator::add_point(Solution point, bool direction)
{
  assert(point.depth() == d_data.d_stage);
  d_points.push_back(point);
  d_directions.push_back(direction);

  GRBLinExpr beta_x;
  GRBLinExpr tau_theta;

  for (size_t stage = 0; stage != d_beta.size(); ++stage)
    beta_x.addTerms(point.d_x[stage].memptr(),
                    d_beta[stage].data(),
                    d_beta[stage].size());

  tau_theta.addTerms(point.d_theta.data(),
                     d_tau.data(),
                     d_tau.size());

  double rhs = d_tau.size() == point.depth() ? 0.0 : dot(d_data.d_costs, point.d_x.back()) + point.d_theta.back();

  if (direction)
    d_mp.addConstr(-beta_x - tau_theta <= rhs);
  else
    d_mp.addConstr(d_alpha - beta_x - tau_theta <= rhs);

  d_mp.update();
}

Solution Enumerator::point()
{
  vvec xvals;
  xvals.reserve(d_x.size());

  for (vvar &xvars: d_x)
  {
    double *vals = d_sp.get(GRB_DoubleAttr_X, xvars.data(), xvars.size());
    xvals.push_back(arma::vec(vals, xvars.size()));
    delete[] vals;
  }

  double *vals = d_sp.get(GRB_DoubleAttr_X, d_theta.data(), d_theta.size());
  vdouble thetavals(vals, vals + d_theta.size());
  delete[] vals;

  return Solution {xvals, thetavals};
}

Solution Enumerator::direction()
{
  GRBModel lp = d_sp;
  GRBVar *vars = lp.getVars();
    // construct lp relaxation
  int nvars = lp.get(GRB_IntAttr_NumVars);
  vector<char> types(nvars, GRB_CONTINUOUS);
  lp.set(GRB_CharAttr_VType, vars, types.data(), nvars);
    // set params
  lp.set(GRB_IntParam_DualReductions, 0);
  lp.set(GRB_IntParam_InfUnbdInfo, 1);
    // obtain unbounded ray
  lp.optimize();
  assert(lp.get(GRB_IntAttr_Status) == 5);
  double *ray = lp.get(GRB_DoubleAttr_UnbdRay, vars, nvars);

  double *copy = ray;
  vvec x;
  for (auto it = d_x.begin(); it != d_x.end(); ++it)
  {
    x.push_back(arma::vec{copy, it->size()});
    copy += it->size();
  }
  vdouble theta{copy, copy + d_theta.size()};

  delete[] ray;
  delete[] vars;

  return Solution {x, theta};
}

void Enumerator::set_sub(Cut &cut)
{
  for (size_t stage = 0; stage != cut.depth(); ++stage)
    d_sp.set(GRB_DoubleAttr_Obj,
             d_x[stage].data(),
             cut.d_beta[stage].memptr(),
             d_x[stage].size());

  d_sp.set(GRB_DoubleAttr_Obj,
           d_theta.data(),
            cut.d_tau.data(),
            cut.depth());

  d_sp.update();
}




void Enumerator::set_bounds(double M)
{
  d_alpha.set(GRB_DoubleAttr_LB, -M);
  d_alpha.set(GRB_DoubleAttr_UB, M);

  for (vvar &beta : d_beta)
  {
    vector<double> lb (beta.size(), -M);
    vector<double> ub (beta.size(), M);
    d_mp.set(GRB_DoubleAttr_LB, beta.data(), lb.data(), lb.size());
    d_mp.set(GRB_DoubleAttr_UB, beta.data(), ub.data(), ub.size());
  }

  vector<double> ub(d_tau.size(), M);
  d_mp.set(GRB_DoubleAttr_UB, d_tau.data(), ub.data(), ub.size());

  d_mp.update();
}