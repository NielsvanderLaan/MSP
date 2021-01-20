#ifndef MSP_STRUCTS_H
#define MSP_STRUCTS_H

#include <armadillo>
#include <vector>
#include "assert.h"

using namespace std;

typedef vector<arma::vec> vvec;
typedef vector<double> vdouble;

struct Solution
{
  vvec d_x;
  vdouble d_theta;

  int depth()
  {
    return d_theta.size();
  }
};

struct Cut
{
  /*
   * cuts are of the form
   * theta_n >= alpha - beta_n x_n - ... - beta_1 x_1
   *                  - tau_n theta_n - ... - tau_1 theta_1
   * d_beta = [beta_n, ..., beta_1]
   * d_tau  = [tau_n, ..., tau_1]
   * note that theta_n features with coefficient 1 + tau_n
   */

  double d_alpha;
  vvec d_beta;
  vdouble d_tau;

  Cut(double a, vvec b, vdouble t)
  :
  d_alpha(a),
  d_beta(b),
  d_tau(t)
  {}

  Cut (vector<int> nvars)
  :
  d_alpha(0)
  {
    for (int n : nvars)
    {
      d_beta.push_back(arma::vec(n, arma::fill::zeros));
      d_tau.push_back(0);
    }
  }

  Cut operator*(double scale)
  {
    vvec beta(this->d_beta);
    for_each(beta.begin(), beta.end(), [scale](arma::vec &vec){ vec *= scale; });

    vdouble tau(this->d_tau);
    for_each(tau.begin(), tau.end(), [scale](double &val){ val *= scale; });

    return Cut{ this->d_alpha * scale, beta, tau};
  }

  Cut& operator+=(const Cut &right)
  {
    this->d_alpha += right.d_alpha;

    transform(this->d_beta.begin(), this->d_beta.end(), right.d_beta.begin(), this->d_beta.begin(), plus<arma::vec>());
    transform(this->d_tau.begin(), this->d_tau.end(), right.d_tau.begin(), this->d_tau.begin(), plus<double>());

    return *this;
  }

  int depth()
  {
    return d_tau.size();
  }
};

double compute_lhs(Cut &cut, Solution &sol);
double compute_rhs(Cut &cut, Solution &sol);
double scaled_rhs(Cut &cut, Solution &sol);

#endif //MSP_STRUCTS_H
