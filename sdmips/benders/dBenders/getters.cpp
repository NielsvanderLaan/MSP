#include "dBenders.h"

Master &dBenders::get_master(int stage, int node)
{
  return d_nodes[stage][node].first;
}

v_enum &dBenders::get_enums(int stage, int node)
{
  return *d_nodes[stage][node].second;
}
