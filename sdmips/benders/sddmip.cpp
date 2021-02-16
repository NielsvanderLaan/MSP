#include "benders.h"

void Benders::sddmip(bool affine)
{
  size_t max_iter = 100;
  for (int iter = 0; iter != max_iter; ++iter)                    // TODO: stopping criterion
  {
    vector<vpath> paths = sample(30);
    vector<vsol> sols = forward(paths, false);

    print_root();

    backward(sols, paths, affine);
  }
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
      sols.push_back(mp.forward());
    }
    ret.emplace_back(move(sols));
  }
  return ret;
}

void Benders::backward(vector<vsol> const &sols, vector<vpath> const &paths, bool affine)
{
  if (d_depth == 0)
    return shared_backward(sols, affine);

  for (int stage = d_data.nstages() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(paths.size());

    for (size_t idx = 0; idx != paths.size(); ++idx)
      cuts.emplace_back(scaled_cut(stage,
                                   master_idx(stage, paths[idx]),
                                   sols[idx][stage],
                                   affine));
        // do not combine the for loops: cuts should be added in bulk
    int count = 0;
    for (size_t idx = 0; idx != cuts.size(); ++idx)
    {
      if (cuts[idx].is_proper(sols[idx][stage]))
      {
        add_cut(cuts[idx], stage, paths[idx]);
        ++count;
      }
    }
    cout << count << " out of " << cuts.size() << " added\n";
  }
  // root node: solutions are identical, just add the cut once
  Cut cut = scaled_cut(0, 0, sols[0][0], affine);
  add_cut_to_root(cut, sols[0][0]);
}

void Benders::shared_backward(vector<vsol> const &sols, bool affine)
{
  for (int stage = d_data.nstages() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(sols.size());
    for (size_t idx = 0; idx != sols.size(); ++idx)
      cuts.emplace_back(shared_scaled_cut(stage,
                                          sols[idx][stage],
                                          affine));

    int count = 0;
    for (size_t idx = 0; idx != cuts.size(); ++idx)
    {
      if (cuts[idx].is_proper(sols[idx][stage]))
      {
        add_shared_cut(cuts[idx], stage);
        ++count;
      }
    }
    cout << count << " out of " << cuts.size() << " added\n";
  }

  Cut cut = shared_scaled_cut(0, sols[0][0], affine);
  add_cut_to_root(cut, sols[0][0]);
}

void Benders::add_cut_to_root(Cut &cut, const Solution &sol)
{
  if (cut.is_proper(sol))
    add_shared_cut(cut, 0);
}

void Benders::print_root()
{
  Master &root = get_master(0, 0);
  root.solve_mip();
  cout << root.mip_obj() << endl;
}
