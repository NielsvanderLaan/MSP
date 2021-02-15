#include "benders.h"

Benders::Benders(GRBEnv &env, Stagewise &data, int depth)
:
d_env(env),
d_data(data),
d_depth(depth),
d_engine(891)                 // random_device{}()
{}