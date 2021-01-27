cmake_minimum_required(VERSION 3.16)
project(MSP)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

find_package(Gurobi REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS -g)

include_directories("C:/armadillo-10.2.0/include")
include_directories("C:/xianyi-OpenBLAS-d2b11c4")
link_directories("C:/gurobi903/win64/lib")

add_executable(MSP main.cpp tree/tree.h tree/add_node.cpp  tree/nodedata/to_model.cpp tree/decomposition/master/master.h tree/decomposition/master/master.cpp tree/decomposition/structs/structs.h tree/decomposition/master/add_cut.cpp tree/decomposition/master/update.cpp tree/decom.cpp tree/decomposition/master/opt_cut.cpp tree/nodedata/nvars.cpp tree/nodedata/ncons.cpp tree/auxiliary/path_to.cpp tree/decomposition/master/master_src.cpp tree/sddp_cut.cpp data/instances.h data/ssv.cpp tree/decomposition/enumerator/enumerator.h tree/decomposition/enumerator/enumerator.cpp tree/auxiliary/nvars.cpp tree/decomposition/structs/src_structs.cpp tree/decomposition/enumerator/add_cut.cpp tree/add_cut.cpp tree/decomposition/enumerator/src_enumerator.cpp tree/decomposition/enumerator/generate_cut.cpp tree/auxiliary/init_enums.cpp tree/scaled_cut.cpp tree/fenchel_cut.cpp tree/forward.cpp tree/auxiliary/is_leaf.cpp tree/solve_master.cpp tree/backward.cpp tree/decomposition/cutfamily/cutfamily.h tree/fp_iteration.cpp tree/lsde.cpp tree/auxiliary/add_to_lsde.cpp tree/solve.cpp tree/cpt_scaled_cut.cpp data/control_1D.cpp)


target_link_libraries(MSP armadillo ${GUROBI_LIBRARIES})
