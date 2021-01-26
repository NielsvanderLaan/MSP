#ifndef MSP_CUTFAMILY_H
#define MSP_CUTFAMILY_H

#include "../structs/structs.h"

class CutFamily
{
public:
  virtual void set(Solution const &sol);

  virtual void add_cut(Cut const &cut);
  virtual Cut compute_cut(double rho, double tol);
  virtual double crho();
};

#endif //MSP_CUTFAMILY_H
