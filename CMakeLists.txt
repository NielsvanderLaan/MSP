cmake_minimum_required(VERSION 3.17)
project(MSP)

set(CMAKE_CXX_STANDARD 14)

include_directories("C:/gurobi903/win64/include")
include_directories("C:/armadillo-10.2.0/include")
include_directories("C:/xianyi-OpenBLAS-d2b11c4")
link_directories("C:/gurobi903/win64/lib")

add_executable(MSP main.cpp tree/tree.h tree/add_node.cpp  tree/nodedata/to_model.cpp lsde/lsde.h lsde/lsde.cpp lsde/add.cpp tree/decomposition/master/master.h tree/decomposition/master/master.cpp tree/decomposition/structs/structs.h tree/decomposition/master/add_cut.cpp tree/decomposition/master/update.cpp tree/decomposition/master/optimize.cpp tree/init_decom.cpp tree/decomposition/master/compute_cut.cpp tree/nodedata/nvars.cpp tree/nodedata/ncons.cpp tree/auxiliary/path_to.cpp tree/decomposition/master/getters.cpp tree/lp_cut.cpp data/instances.h data/ssv.cpp tree/decomposition/enumerator/enumerator.h tree/decomposition/enumerator/enumerator.cpp tree/auxiliary/nvars.cpp tree/decomposition/structs/src.cpp tree/decomposition/enumerator/add_cut.cpp tree/add_cut.cpp tree/decomposition/enumerator/add_cut_to_sp.cpp tree/decomposition/enumerator/add_cut_to_mp.cpp)

target_link_libraries (MSP gurobi_c++)
target_link_libraries (MSP gurobi90)