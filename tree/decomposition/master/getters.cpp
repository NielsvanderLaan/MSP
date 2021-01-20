#include "master.h"

arma::vec Master::xvals()
{
  double *xvals = d_model.get(GRB_DoubleAttr_X, d_xvars.data(), d_xvars.size());
  arma::vec ret(xvals, d_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::theta()
{
  return d_theta.get(GRB_DoubleAttr_X);
}

arma::vec Master::multipliers()
{
  GRBConstr *cons = d_model.getConstrs();
  int ncons = d_data.ncons();

  double *pi = d_model.get(GRB_DoubleAttr_Pi, cons, ncons);
  arma::vec ret(pi, ncons);

  delete[] cons;
  delete[] pi;
  return ret;
}

double Master::obj()
{
  return d_model.get(GRB_DoubleAttr_ObjBound);
}