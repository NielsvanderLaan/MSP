#ifndef MSP_DBENDERS_H
#define MSP_DBENDERS_H

#include "../benders.h"

typedef pair<Master, shared_ptr<v_enum>> node;
typedef vector<node> vnode;

class dBenders : public Benders
{
public:
  vector<vnode> d_nodes;

  dBenders(GRBEnv &env, Stagewise &data, int depth);

  void add_cut(Cut &cut, int stage, int node) override;
  void add_shared_cut(Cut &cut, int stage) override;

  Master &get_master(int stage, int node) override;
  v_enum &get_enums(int stage, int node) override;

};

#endif //MSP_DBENDERS_H
