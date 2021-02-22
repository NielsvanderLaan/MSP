#ifndef MSP_SPBENDERS_H
#define MSP_SPBENDERS_H

#include <memory>
#include "../benders.h"

typedef vector<outer_apx> stage_apx;

class spBenders: public Benders
{
public:
  vector<stage_apx> d_nodal_apx;
  vector<outer_apx> d_stage_apx;

  unique_ptr<Master> d_master;
  v_enum d_gens;

  Master d_root;

  int d_link_depth;

  spBenders(GRBEnv &env, Stagewise &data, int depth, int link_depth);
  spBenders(GRBEnv &env, Stagewise &data, int depth);

  outer_apx const &apx(int stage, int node);

  void add_cut(Cut &cut, int stage, int node) override;
  void add_shared_cut(Cut &cut, int stage) override;

  Master &get_master(int stage, int node) override;
  v_enum &get_enums(int stage, int node) override;

  vector<outer_apx> export_cuts() override;
};

#endif //MSP_SPBENDERS_H
