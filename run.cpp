#include "main.h"

void run(GRBEnv &env,
         Stagewise &problem,
         vector<Family> const &ws_types,
         vector<Family> const &types,
         vector<int> const &iter_limits,
         int depth,
         bool sparse,
         int samples)
{
  if (not ws_types.empty())
  {
    cout << "Warm start: stagewise outer approximations.\nDense implementation.\n";
    cout << "Number of samples: " << samples << '\n';
  }
  dBenders warm_start(env, problem, 0);
  int count = 0;
  for (Family const &type : ws_types)
  {
    cout << to_string(type) << ", " << iter_limits[count] << " iterations.\n";
    warm_start.decom(type, iter_limits[count], false, samples);
    ++count;
  }


  unique_ptr<Benders> benders;
  if (sparse)
    benders = make_unique<spBenders>(env, problem, depth);
  else
    benders = make_unique<dBenders>(env, problem, depth);
  benders->import_cuts(warm_start.export_cuts());

  cout << "SDDMIP, depth: " << depth << ".\n";
  cout << (sparse ? "Sparse" : "Dense")   << " implementation.\n";
  for (Family const &type : types)
  {
    cout << to_string(type) << ", " << iter_limits[count] << " iterations\n";
    benders->decom(type, iter_limits[count], false, samples);
    ++count;
  }

}

void print(int argc, char *argv[])
{
  for (int arg = 1; arg != argc; ++arg)
    cout << argv[arg] << ' ';
  cout << '\n';
}

Args parse(int argc, char *argv[])
{
  string problem = find("PROBLEM=", argc, argv);
  int nstages = stoi(find("STAGES=", argc, argv));
  int n_outcomes = stoi(find("OUTCOMES=", argc, argv));
  int depth = stoi(find("DEPTH=", argc, argv));

  int samples = 30;
  string samples_string = find("SAMPLES=", argc, argv);
  if (valid(samples_string))
    samples = stoi(samples_string);

  bool sparse = valid(find("SPARSE", argc, argv));

  vector<Family> ws_types;
  for (string const &ws : split(find("WSCUTS=", argc, argv), ","))
    ws_types.push_back(to_type(ws));

  vector<Family> types;
  for (string const &ws : split(find("CUTS=", argc, argv), ","))
    types.push_back(to_type(ws));

  vector<int> iter_limits;
  vector<string> iters = split(find("MAX_ITER=", argc, argv), ",");

  if (iters.empty())
    iter_limits = vector<int>(ws_types.size() + types.size(), 25);
  for (string const &limit : iters)
    iter_limits.push_back(stoi(limit));
  cout << "size: " << ws_types.size() << '\n';

  assert(iter_limits.size() == ws_types.size() + types.size());

  return Args{problem, nstages, n_outcomes, ws_types, types, depth, sparse, samples, iter_limits};
}

Stagewise get_problem(Args const &args)
{
  cout << args.problem << ". stages: " << args.nstages << ". outcomes: " << args.n_outcomes << '\n';
  if (args.problem == "CTRL")
    return ctrl_1D(args.nstages, args.n_outcomes);
  if (args.problem == "CLSP")
    return sclsp(args.nstages, args.n_outcomes);

  cout << "unknown problem type\n";
  return Stagewise{};
}




string find(string const &flag, int argc, char *argv[])
{
  for (size_t idx = 1; idx != argc; ++idx)
  {
    string arg(argv[idx]);
    if (arg.find(flag) == 0)
      return arg.substr(flag.size(), arg.size() - flag.size());
  }
  return "NA";
}

Family to_type(string const &type)
{
  if (type == "SDDP")
    return SDDP;
  if (type == "LR")
    return LR;
  if (type == "SC")
   return SC;
  return DEFAULT;
}

vector<string> split(string const &line, string const &delim)
{
  if (not valid(line))
    return vector<string>{};

  vector<string> ret;
  auto start = 0U;
  auto end = line.find(delim, start);
  do
  {
    end = line.find(delim, start);
    ret.push_back(line.substr(start, end - start));
    start = end + delim.length();
  } while (end != std::string::npos);
  return ret;
}

bool valid(string const &test)
{
  return test != "NA";
}




