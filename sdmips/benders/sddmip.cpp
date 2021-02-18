#include "benders.h"

void Benders::decom(Family type, size_t max_iter, bool lp, size_t nsamples)
{
  for (int iter = 0; iter != max_iter; ++iter)                    // TODO: stopping criterion, print more info
  {
    vector<vpath> paths = sample(nsamples);
    vector<vsol> sols = forward(paths, lp);

    print_root();
    backward(type, sols, paths);
  }

  print_root();
}


vector<vsol> Benders::forward(vector<vpath> const &paths, bool lp)
{
  vector<vsol> ret;
  ret.reserve(paths.size());

  for (auto const &path : paths)
  {
    vsol sols;
    sols.reserve(d_data.nstages() - 1);
    for (int stage = 0; stage != d_data.nstages() - 1; ++stage)
    {
      Master &mp = get_master(stage, master_idx(stage, path));
      if (stage > 0)
        mp.update(sols.back());

      lp ? mp.solve_lp() : mp.solve_mip();
      sols.emplace_back(mp.forward());
    }
    ret.emplace_back(move(sols));
  }
  return ret;
}

void Benders::backward(Family type, vector<vsol> const &sols, vector<vpath> const &paths)
{
  bool shared = (type == SDDP) or (d_depth == 0);

  for (int stage = d_data.nstages() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(paths.size());

    for (size_t idx = 0; idx != paths.size(); ++idx)
      cuts.emplace_back(compute_cut(type,
                                    sols[idx][stage],
                                    stage,
                                    shared ? 0 : master_idx(stage, paths[idx])));

        // do not combine the for loops: cut have to be added in bulk
    for (size_t idx = 0; idx != cuts.size(); ++idx)
      add_cut(cuts[idx],
              sols[idx][stage],
              shared,
              stage,
              shared ? -1 : master_idx(stage, paths[idx]));
  }
      // root node: sols are identical
  Cut cut = compute_cut(type, sols[0][0], 0, 0);
  add_cut(cut, sols[0][0], shared, 0, 0);
}

void Benders::add_cut(Cut &cut, Solution const &sol, bool shared, int stage, int node)
{
  if (not shared) assert(node != -1);

  if (cut.is_proper(sol))
    shared ? add_shared_cut(cut, stage) : add_cut(cut, stage, node);
}

void Benders::print_root()
{
  Master &root = get_master(0, 0);
  root.solve_mip();
  cout << root.mip_obj() << endl;
}
