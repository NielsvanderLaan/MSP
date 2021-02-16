#include "benders.h"

void Benders::sddmip(bool affine, size_t nsamples)
{
  size_t max_iter = 100;
  for (int iter = 0; iter != max_iter; ++iter)                    // TODO: stopping criterion
  {
    vector<vpath> paths = sample(nsamples);
    vector<vsol> sols = forward(paths, false);

    print_root();
    backward(sols, paths, affine);
  }

  print_root();
}


void Benders::sddp(size_t nsamples)
{
  size_t max_iter = 5;

  for (int iter = 0; iter != max_iter; ++iter)
  {
    vector<vpath> paths = sample(nsamples);
    vector<vsol> sols = forward(paths, false);

    print_root();
    sddp_backward(sols);
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

    for (size_t idx = 0; idx != cuts.size(); ++idx)
      add_cut(cuts[idx], stage,  sols[idx][stage], paths[idx]);

  }

  Cut cut = scaled_cut(0, 0, sols[0][0], affine);
  add_cut(cut, 0, sols[0][0]);
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

    for (size_t idx = 0; idx != cuts.size(); ++idx)
    {
      if (cuts[idx].is_proper(sols[idx][stage]))
        add_shared_cut(cuts[idx], stage);
    }
  }

  Cut cut = shared_scaled_cut(0, sols[0][0], affine);
  add_cut(cut, 0, sols[0][0]);
}

void Benders::sddp_backward(vector<vsol> const &sols)
{
  for (int stage = d_data.nstages() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(sols.size());
    for (size_t idx = 0; idx != sols.size(); ++idx)
      cuts.emplace_back(sddp_cut(stage, sols[idx][stage]));

    for (size_t idx = 0; idx != cuts.size(); ++idx)
    {
      if (cuts[idx].is_proper(sols[idx][stage]))
        add_shared_cut(cuts[idx], stage);
    }
  }

  Cut cut = sddp_cut(0, sols[0][0]);
  add_cut(cut, 0, sols[0][0]);
}

void Benders::add_cut(Cut &cut, int stage,  Solution const &sol, vpath const &path)
{
  if (cut.is_proper(sol))
    add_cut(cut, stage, path);
}

void Benders::print_root()
{
  Master &root = get_master(0, 0);
  root.solve_mip();
  cout << root.mip_obj() << endl;
}
