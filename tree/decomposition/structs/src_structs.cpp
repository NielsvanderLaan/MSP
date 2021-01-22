#include "structs.h"

int Solution::depth() const
{
  assert(d_x.size() == d_theta.size());
  return d_x.size();
}

void Solution::extend(arma::vec const &x_n, double theta_n)
{
  d_x.push_back(x_n);
  d_theta.push_back(theta_n);
}

Cut::Cut(double a, vvec b, vdouble t)
:
d_alpha(a),
d_beta(b),
d_tau(t)
{}

Cut::Cut (vector<int> nvars)
:
d_alpha(0)
{
  for (int n : nvars)
  {
    d_beta.push_back(arma::vec(n, arma::fill::zeros));
    d_tau.push_back(0);
  }
}

int Cut::depth() const
{
  assert(d_tau.size()== d_beta.size());
  return d_tau.size();
}

Cut operator*(double scale, const Cut &other)
{
  return other * scale;
}


void Cut::scale()
{
  double kappa = 1 + d_tau.back();
  d_alpha /= kappa;
  for_each(d_beta.begin(), d_beta.end(), [kappa](arma::vec &vec){ vec /= kappa; });
  for_each(d_tau.begin(), d_tau.end(), [kappa](double &val){ val /= kappa; });
  d_tau[depth() - 1] = 0;
}

    // free functions
double compute_lhs(Cut const &cut, Solution const &sol)
{
  int idx = cut.depth() - 1;
  return (1 + cut.d_tau[idx]) * sol.d_theta[idx] + arma::dot(cut.d_beta[idx], sol.d_x[idx]);
}

double compute_rhs(Cut const &cut, Solution const &sol)
{
  double ret = cut.d_alpha;

  for (size_t idx = 0; idx < cut.depth() - 1; ++idx)
  {
    ret -= arma::dot(cut.d_beta[idx], sol.d_x[idx]);
    ret -= cut.d_tau[idx] * sol.d_theta[idx];
  }

  return ret;
}

double scaled_rhs(Cut const &cut, Solution const &sol)
{
  double numerator = cut.d_alpha;

  for (size_t idx = 0; idx < cut.depth(); ++idx)
    numerator -= arma::dot(cut.d_beta[idx], sol.d_x[idx]);

  for (size_t idx = 0; idx < cut.depth() - 1; ++idx)
    numerator -= cut.d_tau[idx] * sol.d_theta[idx];

  return numerator / (1 + cut.d_tau.back());
}
















