#include "enumerator.h"

void Enumerator::solve_mp()
{
  d_mp.optimize();
  assert(d_mp.get(GRB_IntAttr_Status) == 2);
}

void Enumerator::solve_sp()
{
  d_sp.optimize();
  assert(d_sp.get(GRB_IntAttr_Status) == 2);
}

double Enumerator::sub_val()
{
  return d_sp.get(GRB_DoubleAttr_ObjVal);
}

double Enumerator::sub_bound()
{
  return d_sp.get(GRB_DoubleAttr_ObjBound);
}

double Enumerator::alpha()
{
  return d_alpha.get(GRB_DoubleAttr_X);
}

double Enumerator::crho()
{
  return d_mp.get(GRB_DoubleAttr_ObjVal) + alpha() - sub_bound();
}

void Enumerator::set_rho(double rho)
{
  d_tau.back().set(GRB_DoubleAttr_Obj, rho);
}

void Enumerator::set_mp(Solution const &sol)
{
  assert(sol.depth() == d_data.d_stage - 1);
  // mind that sol = ([x_a(n)], [theta_a(n)]) (from outside)

  for (size_t stage = 0; stage != d_beta.size(); ++stage)
    d_mp.set(GRB_DoubleAttr_Obj,
             d_beta[stage].data(),
             sol.d_x[stage].memptr(),
             d_beta[stage].size());

  d_mp.set(GRB_DoubleAttr_Obj,
           d_tau.data(),
           sol.d_theta.data(),
           d_tau.size());

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

  return Cut {d_alpha.get(GRB_DoubleAttr_X),
              beta_vals,
              tau_vals};
}

void Enumerator::add_point(Solution point)
{
  //  alpha - beta[x_a(n)] - tau[theta_a(n)] <= c_n x_n + theta_n
  assert(point.depth() == d_data.d_stage);

  d_points.push_back(point);

  GRBLinExpr beta_x;
  GRBLinExpr tau_theta;

  for (size_t stage = 0; stage != d_beta.size(); ++stage)
    beta_x.addTerms(point.d_x[stage].memptr(),
                    d_beta[stage].data(),
                    d_beta[stage].size());

  tau_theta.addTerms(point.d_theta.data(),
                     d_tau.data(),
                     d_tau.size());

  d_mp.addConstr(d_alpha - beta_x - tau_theta <= dot(d_data.d_costs, point.d_x.back()) + point.d_theta.back());
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

void Enumerator::set_sub(Cut &cut)
{
  assert(cut.depth() == d_data.d_stage - 1);

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