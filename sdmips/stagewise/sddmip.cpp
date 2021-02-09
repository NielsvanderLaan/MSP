#include "stagewise.h"

void Stagewise::sddmip(bool affine)
{
  size_t max_iter = 100;
  for (int iter = 0; iter != max_iter; ++iter)                    // TODO: stopping criterion
  {
    vector<vpath> paths = enumerate_paths();      //sample(5);
    vector<vsol> sols = forward(paths, affine, false);
    cout << get_master(0,0).obj() << endl;

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
      int child = master_idx(stage + 1, path);

      solve(stage, master, affine, lp, true);
      Solution forward = solution(stage, master);
      sols.push_back(forward);

      get_master(stage + 1, child).update(forward);
    }
    ret.emplace_back(move(sols));
  }

  return ret;
}

void Stagewise::backward(vector<vsol> const &sols, vector<vpath> const &paths, bool affine)
{
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

void Stagewise::solve(int stage, int node, bool affine, bool lp, bool force)
{
  Master &mp = get_master(stage, node);
  mp.solve_lp();

  if (lp)
    return;

  while (not mp.integer())
  {
    Cut cutting_plane = fenchel_cut(stage, node, affine);
    if (not add_cp(cutting_plane, stage, node))    // mip could not be solved using cutting planes
    {
      if (force)
        mp.solve_mip();                              // use Gurobi
      return;
    }

    mp.solve_lp();
  }
}


