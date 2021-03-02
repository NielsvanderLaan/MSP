#ifndef MSP_DBENDERS_H
#define MSP_DBENDERS_H

#include "../benders.h"

typedef pair<Master, shared_ptr<v_enum>> node;
typedef vector<node> vnode;
typedef vector<Gomory> vgom;

class dBenders : public Benders
{
public:
  vector<vnode> d_nodes;
  vector<vgom> d_gomory;

  dBenders(GRBEnv &env, Stagewise &data, int depth, int link_depth);
  dBenders(GRBEnv &env, Stagewise &data, int depth);

  void init_gomory(GRBEnv &env);
  void add_cut(Cut &cut, int stage, int node) override;
  void add_shared_cut(Cut &cut, int stage) override;

  Master &get_master(int stage, int node) override;
  v_enum &get_enums(int stage, int node) override;
  Gomory &get_gomory(int stage, int out) override;

  vector<outer_apx> export_cuts() override;
};

#endif //MSP_DBENDERS_H
