#ifndef MSP_STRUCTS_H
#define MSP_STRUCTS_H

#include <armadillo>
#include <vector>
#include "assert.h"
#include "gurobi_c++.h"

using namespace std;

typedef vector<arma::vec> vvec;
typedef vector<double> vdouble;

struct Solution
{
  vvec d_x;         // (x_1,...,x_a(n), x_n)
  vdouble d_theta;  // (theta_1,...,theta_a(n), theta_n)

  int depth() const;
  void extend(arma::vec const &x_n, double theta_n);

  void print();
};

struct Cut
{
  /*
   * cuts are of the form
   * theta_n >= alpha - beta_n x_n - ... - beta_1 x_1
   *                  - tau_n theta_n - ... - tau_1 theta_1
   * d_beta = [beta_1, ..., beta_n]
   * d_tau  = [tau_1, ..., tau_n]
   * note that theta_n features with coefficient 1 + tau_n
   */

  double d_alpha;
  vvec d_beta;
  vdouble d_tau;
  bool d_feas;

  Cut();
  Cut(double a, vvec b, vdouble t);
  Cut (vector<int> const &nvars);


  Cut operator*(double scale) const
  {
    vvec beta(this->d_beta);
    for_each(beta.begin(), beta.end(), [scale](arma::vec &vec){ vec *= scale; });

    vdouble tau(this->d_tau);
    for_each(tau.begin(), tau.end(), [scale](double &val){ val *= scale; });

    return Cut{ this->d_alpha * scale, beta, tau};
  }

  friend Cut operator*(double scale, Cut const &other);

  Cut& operator+=(const Cut &right)
  {
    this->d_alpha += right.d_alpha;

    transform(this->d_beta.begin(), this->d_beta.end(), right.d_beta.begin(), this->d_beta.begin(), plus<arma::vec>());
    transform(this->d_tau.begin(), this->d_tau.end(), right.d_tau.begin(), this->d_tau.begin(), plus<double>());

    return *this;
  }

  int depth() const;
  void scale();
  void print();
};

double compute_lhs(Cut const &cut, Solution const &sol);
double compute_rhs(Cut const &cut, Solution const &sol);
double scaled_rhs(Cut const &cut, Solution const &sol);
bool is_integer(arma::vec x, arma::Col<char> const &types);

#endif //MSP_STRUCTS_H
