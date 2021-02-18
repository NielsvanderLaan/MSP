#include "instances.h"

Stagewise sclsp(size_t nstages, size_t n_outcomes)
{
  int seed = 1234; //random_device{}()
  mt19937 engine(seed);
  lognormal_distribution<double> log_normal(0.0, 1.0);

  int nvars = 12;
  int ncons = 7;
  double I_bar = 1000;
  double K = 175;
  double b_i = 5.0;     // lost sales penalty
  double inf = GRB_INFINITY;
  inf = 1e5;
  double bigM = 1e6;

  arma::vec d1 {40, 35, 0};
  arma::vec d2 {45, 35, 60};
  arma::vec d3 {65, 55, 45};
  arma::vec d4 {35, 35, 80};
  arma::mat Demand = join_rows(d1, d2, d3, d4);
  arma::vec mean_demand = mean(Demand, 1);


  arma::vec setup_time {15, 10, 15};
  arma::vec prod_time {1, 1 ,1};

    // variables are ordered as (S, P, I, L)
  arma::vec costs {60,  120,  80,      // setup
                   0,   0 ,   0,       // unit production
                   2,   3,    1,       // holding costs
                   b_i, b_i,  b_i};    // lateness costs

 arma::vec lb(nvars, arma::fill::zeros);
 arma::vec ub {1,     1,     1,          // setup (binary)
               inf,   inf,   inf,        // production
               I_bar, I_bar, I_bar,      // inventory capacity
               inf,   inf,   inf};       // lost sales


  arma::sp_mat Amat(nvars, ncons);
  arma::sp_mat Bmat(nvars, ncons);

  for (int i = 0; i != Demand.n_rows; ++i)
  {
    Amat(i + 3, i) =  1;      // P_{t, i}
    Amat(i + 6, i) = -1;      // -I_{t, i}
    Amat(i + 9, i) =  1;      // +L_{t, i}
    Bmat(i + 6, i) =  1;      // +I_{t - 1, i}

    Amat(i,            i + 3) = -bigM;    // -M*S_i
    Amat(i + 3, i + 3) = 1;        // +P_i
  }
  Amat.submat(0, 6, 5, 6) = arma::join_cols(setup_time, prod_time);

  arma::vec rhs = arma::join_cols(Demand.col(0),
                                  arma::vec(3, arma::fill::zeros),
                                  arma::vec{K});
  vector<char> types(nvars, GRB_CONTINUOUS);
  fill_n(types.begin(), 3, GRB_INTEGER);

  vector<char> senses(ncons, GRB_LESS_EQUAL);
  fill_n(senses.begin(), 3, GRB_EQUAL);

  NodeData root{1,
                1.0,
                0.0,
                costs,
                lb,
                ub,
                Amat,
                sp_mat(),
                rhs,
                types,
                senses,
                {}};

  Stagewise sw;
  sw.add_node(root);

  NodeData sub = root;
  sub.d_Bmat = Bmat;
  sub.d_prob = 1.0 / n_outcomes;

  for (size_t stage = 1; stage != nstages; ++stage)
  {
    sub.d_stage = stage + 1;
    for (int out = 0; out != n_outcomes; ++out)
    {
      arma::vec base_demand = stage < Demand.n_cols ? Demand.col(stage) : mean_demand;

      for (int i = 0; i != base_demand.n_rows; ++i)    // product type
        sub.d_rhs[i] = base_demand(i) * log_normal(engine);

      sw.add_node(sub);
    }
  }

  return sw;
}