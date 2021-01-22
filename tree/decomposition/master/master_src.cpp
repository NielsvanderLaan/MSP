#include "master.h"

void Master::solve_lp()
{
  d_lp.optimize();
}

void Master::solve_mip()
{
  d_mip.optimize();
}

arma::vec Master::mip_xvals()
{
  double *xvals = d_mip.get(GRB_DoubleAttr_X, d_mip_xvars.data(), d_mip_xvars.size());
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
  double *xvals = d_lp.get(GRB_DoubleAttr_X, d_lp_xvars.data(), d_lp_xvars.size());
  arma::vec ret(xvals, d_lp_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::lp_theta()
{
  return d_lp_theta.get(GRB_DoubleAttr_X);
}

arma::vec Master::multipliers()
{
  GRBConstr *cons = d_lp.getConstrs();
  int ncons = d_data.ncons();

  double *pi = d_lp.get(GRB_DoubleAttr_Pi, cons, ncons);
  arma::vec ret(pi, ncons);

  delete[] cons;
  delete[] pi;
  return ret;
}

double Master::obj()
{
  return d_lp.get(GRB_DoubleAttr_ObjBound);
}