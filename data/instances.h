#ifndef MSP_INSTANCES_H
#define MSP_INSTANCES_H

#include <random>
#include "../sdmips/tree/tree.h"
#include "../sdmips/stagewise/stagewise.h"

using namespace std;
using namespace arma;

Tree ssv();
Tree control_1D();
Stagewise ctrl_1D();

#endif //MSP_INSTANCES_H
