#include "master.h"

void Master::solve_lp()
{
  d_lp->optimize();
  int status = d_lp->get(GRB_IntAttr_Status);

  assert(d_lp->get(GRB_IntAttr_Status) == 2);
  d_x_n = lp_xvals();
  d_theta_n = lp_theta();
}

void Master::solve_mip()
{
  d_mip->optimize();
  assert(d_mip->get(GRB_IntAttr_Status) == 2);
  d_x_n = mip_xvals();
  d_theta_n = mip_theta();
}

arma::vec Master::mip_xvals()
{
  double *xvals = d_mip->get(GRB_DoubleAttr_X, d_mip_xvars.data(), d_mip_xvars.size());
  arma::vec ret(xvals, d_mip_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::mip_theta()
{
  return d_mip_theta.get(GRB_DoubleAttr_X);
}

arma::vec Master::lp_xvals()
{
  double *xvals = d_lp->get(GRB_DoubleAttr_X, d_lp_xvars.data(), d_lp_xvars.size());
  arma::vec ret(xvals, d_lp_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::lp_theta()
{
  return d_lp_theta.get(GRB_DoubleAttr_X);
}

Solution Master::forward()
{
  Solution sol = d_state;
  sol.extend(d_x_n, d_theta_n);
  return sol;
}

double Master::theta_n()
{
  return d_theta_n;
}

bool Master::integer()
{
  return is_integer(lp_xvals(), d_data.d_types);
}

arma::vec Master::multipliers()
{
  GRBConstr *cons = d_lp->getConstrs();
  int ncons = d_data.ncons() + d_cuts.size();

  double *pi = d_lp->get(GRB_DoubleAttr_Pi, cons, ncons);
  arma::vec ret(pi, ncons);

  delete[] cons;
  delete[] pi;
  return ret;
}

double Master::obj()
{
  return d_lp->get(GRB_DoubleAttr_ObjBound);
}

void Master::set_rho(double rho)
{
  double diff = rho - d_state.d_theta.back() ;
  d_state.d_theta.back() = rho;
  GRBConstr *cons = d_lp->getConstrs();
  double *rhs = d_lp->get(GRB_DoubleAttr_RHS, cons + d_data.ncons(), d_cuts.size());

  for (size_t cut = 0; cut != d_cuts.size(); ++cut)
  {
    vdouble const &tau = d_cuts[cut].d_tau;
    rhs[cut] -= diff * tau[tau.size() - 2];
  }

  d_lp->set(GRB_DoubleAttr_RHS, cons + d_data.ncons(), rhs, d_cuts.size());
  d_lp->update();

  delete[] cons;
  delete[] rhs;
}