#include "enumerator.h"

void Enumerator::optimize_mp()
{
  d_mp->optimize();

  if (mp_status() != 2)
    throw mp_exception{};
}

void Enumerator::solve_sp()
{
  d_sp->optimize();

  int status = sp_status();
  if (status == GRB_OPTIMAL or status == GRB_SUBOPTIMAL)
    return;

  if (status == GRB_TIME_LIMIT)
  {
    try
    {
      sub_val();
    } catch (GRBException &e)
    {
      throw sp_exception{};
    }
    return;
  }

  throw sp_exception{};
}

int Enumerator::sp_status() const
{
  return d_sp->get(GRB_IntAttr_Status);
}

int Enumerator::mp_status() const
{
  return d_mp->get(GRB_IntAttr_Status);
}

double Enumerator::mp_violation() const
{
  return d_mp->get(GRB_DoubleAttr_ConstrVio) + d_mp->get(GRB_DoubleAttr_ConstrResidual);
}

double Enumerator::sub_val() const
{
  return d_sp->get(GRB_DoubleAttr_ObjVal);
}

double Enumerator::sub_bound() const
{
  if (sp_status() != 2)
    return -GRB_INFINITY;

  return d_sp->get(GRB_DoubleAttr_ObjBound);
}

double Enumerator::alpha() const
{
  return d_alpha.get(GRB_DoubleAttr_X);
}

double Enumerator::crho() const
{
  if (mp_status() != 2)
    return GRB_INFINITY;

  return -d_mp->get(GRB_DoubleAttr_ObjVal) + alpha() - sub_bound();
}

void Enumerator::set_rho(double rho)
{
  d_mp->chgCoeff(d_objcon, d_tau.back(), rho);
  d_mp->update();
}

void Enumerator::set_mp(Solution const &sol)
{
  for (size_t stage = 0; stage != d_beta.size(); ++stage)
  {
    for (size_t var = 0; var != d_beta[stage].size(); ++var)
      d_mp->chgCoeff(d_objcon, d_beta[stage][var], sol.d_x[stage][var]);
  }

  for (size_t var = 0; var != d_tau.size(); ++var)
    d_mp->chgCoeff(d_objcon, d_tau[var], sol.d_theta[var]);

  d_mp->update();
}

Cut Enumerator::solve_mp(bool affine, double M)
{
  optimize_mp();

  Cut ret = candidate();
  if (ret.abs_max() < M)
    return ret;

  double old_M = d_alpha.get(GRB_DoubleAttr_UB);
  set_bounds(affine, M);
  try
  {
    optimize_mp();
  } catch (mp_exception)
  {
    set_bounds(affine, old_M);
    throw;
  }

  set_bounds(affine, old_M);
  return candidate();
}

Cut Enumerator::candidate()
{
  vvec beta_vals;
  beta_vals.reserve(d_beta.size());

  for (auto &beta : d_beta)
  {
    double *vals = d_mp->get(GRB_DoubleAttr_X, beta.data(), beta.size());
    beta_vals.push_back(arma::vec(vals, beta.size()));
    delete[] vals;
  }

  double *vals = d_mp->get(GRB_DoubleAttr_X, d_tau.data(), d_tau.size());
  vdouble tau_vals(vals, vals + d_tau.size());
  delete[] vals;

  for_each(tau_vals.begin(), tau_vals.end(), [](double &val){ val = max(0.0, val);});

  return Cut {alpha(),
              beta_vals,
              tau_vals};
}

void Enumerator::add_point(Solution point, bool prime)
{
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

  double rhs = d_tau.size() == point.depth() ? 0.0 : dot(d_data.d_costs, point.d_x.back()) + point.d_theta.back();
  d_mp->addConstr(d_alpha - beta_x - tau_theta <= rhs);

  if (prime)
    d_obj.set(GRB_DoubleAttr_UB, rhs);

  d_mp->update();
}

Solution Enumerator::point()
{
  vvec xvals;
  xvals.reserve(d_x.size());

  for (vvar &xvars: d_x)
  {
    double *vals = d_sp->get(GRB_DoubleAttr_X, xvars.data(), xvars.size());
    xvals.push_back(arma::vec(vals, xvars.size()));
    delete[] vals;
  }

  double *vals = d_sp->get(GRB_DoubleAttr_X, d_theta.data(), d_theta.size());
  vdouble thetavals(vals, vals + d_theta.size());
  delete[] vals;

  return Solution {xvals, thetavals};
}

void Enumerator::set_sub(Cut &cut)
{
  for (size_t stage = 0; stage != cut.depth(); ++stage)
    d_sp->set(GRB_DoubleAttr_Obj,
             d_x[stage].data(),
             cut.d_beta[stage].memptr(),
             d_x[stage].size());

  d_sp->set(GRB_DoubleAttr_Obj,
           d_theta.data(),
            cut.d_tau.data(),
            cut.depth());

  d_sp->update();
}

void Enumerator::set_bounds(bool affine, double M)
{
  d_alpha.set(GRB_DoubleAttr_LB, -M);
  d_alpha.set(GRB_DoubleAttr_UB, M);

  for (vvar &beta : d_beta)
  {
    vector<double> lb (beta.size(), -M);
    vector<double> ub (beta.size(), M);
    d_mp->set(GRB_DoubleAttr_LB, beta.data(), lb.data(), lb.size());
    d_mp->set(GRB_DoubleAttr_UB, beta.data(), ub.data(), ub.size());
  }

  set_tau_bounds(affine, M);
}

void Enumerator::set_tau_bounds(bool affine, double M)
{
  vector<double> ub(d_tau.size(), affine ? 0.0 : M);
  d_mp->set(GRB_DoubleAttr_UB, d_tau.data(), ub.data(), ub.size());
}

void Enumerator::set_mp(bool tight)
{
  d_mp->set(GRB_IntParam_ScaleFlag, tight ? 0 : -1);
  d_mp->set(GRB_IntParam_NumericFocus, tight ? 3: 0);
  d_mp->set(GRB_IntParam_Presolve, tight ? 0 : -1);
  d_mp->reset();
}

void Enumerator::set_sp(bool tight)
{
  d_sp->set(GRB_IntParam_Method, tight ? 0 : -1);
  d_sp->set(GRB_IntParam_NumericFocus, tight ? 3: 0);
  d_sp->reset();
}

void Enumerator::clear()
{
  assert(d_mp->get(GRB_IntAttr_NumConstrs) == d_points.size() + 1);

  GRBConstr *cons = d_mp->getConstrs();
  for (size_t idx = 1; idx != d_mp->get(GRB_IntAttr_NumConstrs); ++idx)
    d_mp->remove(cons[idx]);
  delete[] cons;

  d_points.clear();

  d_mp->update();
}