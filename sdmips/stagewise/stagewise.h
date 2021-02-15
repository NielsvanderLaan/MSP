#ifndef MSP_STAGEWISE_H
#define MSP_STAGEWISE_H

#include <assert.h>
#include "../nodedata/nodedata.h"

typedef vector<NodeData> stage_data;

class Stagewise
{
public:
  vector<stage_data> d_stages;

  void add_node(NodeData const &data);

  GRBModel lsde(GRBEnv &env);

  int nstages() const;
  int outcomes(int stage) const;
  vector<int> nvars(int stage) const;
  vector<double> probs(int stage) const;
};

#endif //MSP_STAGEWISE_H
