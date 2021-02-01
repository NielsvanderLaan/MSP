#include "structs.h"

int Solution::depth() const
{
  assert(d_x.size() == d_theta.size());
  return d_theta.size();
}

double Solution::theta_n() const
{
  return d_theta.back();
}

void Solution::extend(arma::vec const &x_n, double theta_n)
{
  d_x.push_back(x_n);
  d_theta.push_back(theta_n);
}

void Solution::print() const
{
  for (int stage = 0; stage != d_x.size(); ++stage)
  {
    auto &x = d_x[stage];
    cout << "x_" << stage + 1 << ' ';
    for_each(x.begin(), x.end(), [](double val){cout << val << ' ';});
    cout << '\n';
  }
  cout << "theta: ";
  for_each(d_theta.begin(), d_theta.end(),[](double val){cout << val << ' ';});
  cout << '\n';
}

bool is_integer(double val, double precision = 1e-6)
{
  return abs(val - floor(val + 0.5)) <= precision;
}

bool is_integer(arma::vec x, arma::Col<char> const &types)
{
  for (size_t idx = 0; idx != types.size(); ++idx)
  {
    if (types[idx] == GRB_CONTINUOUS)
      continue;
    if (not is_integer(x[idx]))
      return false;
  }

  return true;
}

Cut::Cut()
:
d_feas(false)
{}

Cut::Cut(double a, vvec b, vdouble t)
:
d_alpha(a),
d_beta(move(b)),
d_tau(move(t)),
d_feas(false)
{}

Cut::Cut (vector<int> const &nvars)
:
d_alpha(0),
d_feas(false)
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

void Cut::print() const
{
  cout << "alpha: " << d_alpha;
  for (size_t idx = 0; idx != d_beta.size(); ++idx)
  {
    cout << "\nbeta_" << idx + 1 << ": ";
    for_each(d_beta[idx].begin(), d_beta[idx].end(), [](double val){cout << val << ' ';});
  }
  cout << "\ntau: ";
  for_each(d_tau.begin(), d_tau.end(), [](double val){cout << val << ' ';});
  cout << '\n';
}

void Cut::scale()
{
  double kappa = 1 + d_tau.back();
  d_alpha /= kappa;
  for_each(d_beta.begin(), d_beta.end(), [kappa](arma::vec &vec){ vec /= kappa; });
  for_each(d_tau.begin(), d_tau.end(), [kappa](double &val){ val /= kappa; });
  d_tau.back() = 0;
}

    // free functions
double compute_lhs(Cut const &cut, Solution const &sol)
{
  int idx = cut.depth() - 1;
  return (1 + cut.d_tau.back()) * sol.d_theta[idx] + arma::dot(cut.d_beta.back(), sol.d_x[idx]);
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
















