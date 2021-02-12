#include "stagewise.h"

void Stagewise::sddmip(bool affine)
{
  size_t max_iter = 2;
  for (int iter = 0; iter != max_iter; ++iter)                    // TODO: stopping criterion
  {
    vector<vpath> paths = sample(30);                    // enumerate_paths();
    vector<vsol> sols = forward(paths, affine, false);
    cout << get_master(0,0).mip_obj() << endl;
    backward(sols, paths, affine);
  }
}

vector<vsol> Stagewise::forward(vector<vpath> const &paths, bool affine, bool lp)
{
  vector<vsol> ret;
  ret.reserve(paths.size());

  for (auto const &path : paths)
  {
    vsol sols;
    sols.reserve(d_stages.size() - 1);
    for (int stage = 0; stage != d_stages.size() - 1; ++stage)
    {
      int master = master_idx(stage, path);
      lp ? solve_lp(stage, master) : solve_mip(stage, master);

      Solution forward = solution(stage, master);
      sols.push_back(forward);

      if (stage + 1 == d_stages.size() - 1)     // penultimate stage: no need to update children
        continue;

      int child = master_idx(stage + 1, path);
      get_master(stage + 1, child).update(forward);
    }
    ret.emplace_back(move(sols));
  }
  return ret;
}

void Stagewise::backward(vector<vsol> const &sols, vector<vpath> const &paths, bool affine)
{
  if (d_depth == 0)
    return shared_backward(sols, affine);

  for (int stage = d_stages.size() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(paths.size());

    for (size_t idx = 0; idx != paths.size(); ++idx)
      cuts.emplace_back(scaled_cut(stage,
                                  master_idx(stage, paths[idx]),
                                   sols[idx][stage],
                                   affine));
        // do not combine the for loops: cuts should be added in bulk
    for (size_t idx = 0; idx != paths.size(); ++idx)
      add_cut(cuts[idx], stage, paths[idx]);
  }
        // root node: solutions are identical, just add the cut once
  Cut cut = scaled_cut(0, 0, sols[0][0], affine);
  add_cut(cut, 0, paths[0]);
}

void Stagewise::shared_backward(vector<vsol> const &sols, bool affine)
{
  for (int stage = d_stages.size() - 2; stage != 0; --stage)
  {
    vector<Cut> cuts;
    cuts.reserve(sols.size());
    for (size_t idx = 0; idx != sols.size(); ++idx)
      cuts.emplace_back(shared_scaled_cut(stage,
                                          sols[idx][stage],
                                          affine));
    for (Cut &cut : cuts)
      add_shared_cut(cut, stage);
  }

  Cut cut = shared_scaled_cut(0, sols[0][0], affine);
  add_shared_cut(cut, 0);
}

void Stagewise::solve_lp(int stage, int node)
{
  get_master(stage, node).solve_lp();
}

void Stagewise::solve_mip(int stage, int node)
{
  get_master(stage, node).solve_mip();
}


