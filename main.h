#ifndef MSP_MAIN_H
#define MSP_MAIN_H

#include <iostream>

#include "gurobi_c++.h"

#include "sdmips/tree/tree.h"
#include "sdmips/benders/dbenders/dbenders.h"
#include "sdmips/benders/spbenders/spbenders.h"
#include "data/instances.h"
#include "sdmips/benders/families.h"

using namespace std;
using namespace arma;

void run(GRBEnv &env,
         Stagewise &problem,
         vector<Family> const &ws_types,
         vector<Family> const &types,
         vector<int> const &iter_limits,
         int depth,
         bool sparse,
         int samples);

struct Args
{
  string problem;
  int nstages;
  int n_outcomes;

  vector<Family> ws_types;
  vector<Family> types;

  int depth;
  bool sparse;
  int samples;
  vector<int> iter_limits;
};

Args parse(int argc, char *argv[]);
void print(int argc, char *argv[]);
Stagewise get_problem(Args const &args);

string find(string const &flag, int argc, char *argv[]);
Family to_type(string const &type);
vector<string> split(string const &line, string const &delim);
bool valid(string const &test);








#endif //MSP_MAIN_H
