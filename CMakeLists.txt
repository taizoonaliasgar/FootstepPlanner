cmake_minimum_required(VERSION 3.12)
project(footstep_planner)


set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-O3")
set(USE_QPSWIFT TRUE) 
add_definitions(-DUSE_QPSWIFT=1)


# set(LL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../A1_LL_Exp-Arch_Change)
set(LL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LL_integrator)
set(ROBOT_DIR "A1_New_wFoot")

find_package(raisim CONFIG HINTS "/home/taizoon/raisimEnvT/rInstall/lib/cmake/raisim")
find_package(raisimOgre CONFIG HINTS "/home/taizoon/raisimEnvT/rInstall/share/raisimOgre/cmake")

# find_package(raisim CONFIG HINTS "/home/trec/WorkRaj/raisim_legged/raisim_build/lib/cmake/raisim")
# find_package(OGRE CONFIG REQUIRED PATHS "/home/trec/WorkRaj/raisim_legged/raisim_build/lib/OGRE/cmake")
# find_package(raisimOgre CONFIG HINTS  "/home/trec/WorkRaj/raisim_legged/raisim_build/share/raisimOgre/cmake")  

include_directories(PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include)
find_package(Eigen3 REQUIRED HINTS ${Eigen3_HINT})

#find_package(osqp REQUIRED)

find_package(casadi REQUIRED)
find_package(OpenSSL REQUIRED)
# set casadi include directory
set(CASADI_INCLUDE_DIR /usr/local/include/casadi)

file(GLOB HEADER_FILES_HPP ${CMAKE_SOURCE_DIR}/include/*.hpp)
file(GLOB HEADER_FILES_H ${CMAKE_SOURCE_DIR}/include/*.h)


# find casadi library
find_library(CASADI_LIBRARY
    NAMES casadi
    HINTS ${CASADI_INCLUDE_DIR}/../lib $ENV{CASADI_PREFIX}/lib)
if(CASADI_LIBRARY)
    set(CASADI_LIBRARIES ${CASADI_LIBRARIES} ${CASADI_LIBRARY})
endif()

find_package(Threads REQUIRED)
set(IPOPT_INCLUDE_DIRS /usr/include/coin)

if (Eigen3_FOUND)
    message("Eigen found")
    message("  Eigen location: ${Eigen3_DIR}")
    message("  Eigen Version: ${Eigen3_VERSION}")
endif ()

#HDSRL portion
set(UT_INC ${CMAKE_CURRENT_SOURCE_DIR}/unitree_legged_sdk/include)
set(UT_LINK ${CMAKE_CURRENT_SOURCE_DIR}/unitree_legged_sdk/lib)
set(EXTRA_LIBS -pthread ${UT_LINK}/libunitree_legged_sdk_amd64.so)

include_directories(PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})
#include_directories(${Boost_INCLUDE_DIRS})
set(BOOST_ROOT "/opt/boost_1_68_0")
set(Boost_NO_SYSTEM_PATHS TRUE)
find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})


set(dyn_srcs 
    "${LL_DIR}/dynamics/A1/A1_Dynamics_z.cpp"
    "${LL_DIR}/dynamics/RobotModel.cpp")
    
add_library(dyn_lib ${dyn_srcs})
target_include_directories(dyn_lib PUBLIC
    ${LL_DIR}/dynamics
    ${LL_DIR}/dynamics/A1
    #${LL_DIR}/dynamics/${ROBOT_DIR}
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include) 

set(dyn_srcs_o 
    "${LL_DIR}/dynamics/A1/A1_Dynamics.cpp"
    "${LL_DIR}/dynamics/RobotModel.cpp")
    
add_library(dyn_lib_o ${dyn_srcs_o})
target_include_directories(dyn_lib_o PUBLIC
    ${LL_DIR}/dynamics
    ${LL_DIR}/dynamics/A1
    #${LL_DIR}/dynamics/${ROBOT_DIR}
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)

set(dyn_srcs_full 
    "${LL_DIR}/dynamics/A1/A1_Dynamics_z_full.cpp"
    "${LL_DIR}/dynamics/RobotModel.cpp")
    
add_library(dyn_lib_full ${dyn_srcs_full})
target_include_directories(dyn_lib_full PUBLIC
    ${LL_DIR}/dynamics
    ${LL_DIR}/dynamics/A1
    #${LL_DIR}/dynamics/${ROBOT_DIR}
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include) 
 

# ===== Parameter library ===== #
add_library(params_lib
    "${CMAKE_CURRENT_SOURCE_DIR}/params/Parameters.cpp")
target_include_directories(params_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/params
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)

# ===== Optimization library ===== #
file(GLOB OPTIMIZATION_SRC
    "optimization/iSWIFT/src/Auxilary.c"
    "optimization/iSWIFT/src/Prime.c"
    "optimization/iSWIFT/src/timer.c"
    "optimization/iSWIFT/ldl/src/ldl.c"
    "optimization/iSWIFT/cpp_wrapper/iswift_qp.cpp"
)
add_library(optimization_lib
    ${OPTIMIZATION_SRC}
    )
target_include_directories(optimization_lib PUBLIC
    optimization/iSWIFT/include
    optimization/iSWIFT/ldl/include
    optimization/iSWIFT/cpp_wrapper)

# ===== Control ===== #
add_library(ctrl_lib
    #"${LL_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/VirtualConstraints.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/DataLog.cpp")
target_include_directories(ctrl_lib PUBLIC 
    #${LL_DIR}/include
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
target_link_libraries(ctrl_lib PUBLIC 
    #dyn_lib_full
    dyn_lib
    params_lib
    optimization_lib)

add_library(ctrl_lib_o
    #"${LL_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/VirtualConstraints.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/DataLog.cpp")
target_include_directories(ctrl_lib_o PUBLIC 
    #${LL_DIR}/include
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
target_link_libraries(ctrl_lib_o PUBLIC 
    dyn_lib_o#full
    params_lib
    optimization_lib)

add_library(ctrl_lib_full
    #"${LL_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LowLevelCtrl.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/VirtualConstraints.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/DataLog.cpp")
target_include_directories(ctrl_lib_full PUBLIC 
    #${LL_DIR}/include
    #${LL_DIR}/util_include
    #${LL_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
target_link_libraries(ctrl_lib_full PUBLIC 
    dyn_lib_full
    params_lib
    optimization_lib)

# ===== Planning and Contact ===== #
add_library(plan_lib
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ContactEst.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/MotionPlanner.cpp")
target_include_directories(plan_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    #${LL_DIR}/global_include
    #${LL_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include)

# ===== Upright walk wrapper Library ===== #
# add_library(wrapperwalk
#     "${CMAKE_CURRENT_SOURCE_DIR}/src/LocoWrapperwalk.cpp")
# target_include_directories(wrapperwalk PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}/include)
 
#  target_link_libraries(wrapperwalk PUBLIC
#     ctrl_lib
#     plan_lib)

# add_library(wrapper
#     "${CMAKE_CURRENT_SOURCE_DIR}/src/LocoWrapper.cpp")
# target_include_directories(wrapper PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}/include)
 
#  target_link_libraries(wrapper PUBLIC
#     ctrl_lib_o
#     plan_lib)

# # ===== Combining dynamics for both configurations ===== #
add_library(wrapperwalk2
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LocoWrapperwalk.cpp")
target_include_directories(wrapperwalk2 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
 
 target_link_libraries(wrapperwalk2 PUBLIC
    ctrl_lib_full
    plan_lib)

add_library(wrapper2
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LocoWrapper.cpp")
target_include_directories(wrapper2 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
 
 target_link_libraries(wrapper2 PUBLIC
    ctrl_lib_full
    plan_lib)

# # ===== NMPC Library ===== #
add_library(nmpc_lib
    "src/SRBNMPC.cpp")
target_compile_options(nmpc_lib PRIVATE -std=c++17)

target_include_directories(nmpc_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
target_link_libraries(nmpc_lib PUBLIC 
    params_lib
    -L/usr/local/lib)

add_library(nmpcr_lib
    "src/SRBNMPCR.cpp")
target_compile_options(nmpcr_lib PRIVATE -std=c++17)

target_include_directories(nmpcr_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/util_include
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
target_link_libraries(nmpcr_lib PUBLIC 
    params_lib
    -L/usr/local/lib)

# # ===== Microstrain ===== #
set(MSCL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../MSCL")
# Add the static library
add_library(mscl STATIC IMPORTED)
set_target_properties(mscl PROPERTIES 
    IMPORTED_LOCATION 
    "${MSCL_DIR}/build/MSCL/MSCL_Static/libMSCL.a")
# Include the header files
target_include_directories(mscl INTERFACE 
    "${MSCL_DIR}/MSCL/source")
target_link_libraries(mscl INTERFACE ${Boost_LIBRARIES} OpenSSL::SSL OpenSSL::Crypto)


add_executable(upright_NMPC src/upright_NMPC.cpp)
target_link_libraries(upright_NMPC PUBLIC raisim::raisimOgre raisim::assimp nmpc_lib ${CASADI_LIBRARIES} ipopt)
target_include_directories(upright_NMPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
target_include_directories(upright_NMPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(upright_NMPC PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(upright_NMPC PRIVATE -std=c++17)

# add_executable(upright_NMPCR src/upright_NMPCR.cpp)
# target_link_libraries(upright_NMPCR PUBLIC raisim::raisimOgre raisim::assimp nmpcr_lib ${CASADI_LIBRARIES} ipopt)
# target_include_directories(upright_NMPCR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
# target_include_directories(upright_NMPCR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_compile_definitions(upright_NMPCR PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
# target_compile_options(upright_NMPCR PRIVATE -std=c++17)

add_executable(uprightwalk src/A1_up_walk.cpp)
target_link_libraries(uprightwalk PUBLIC raisim::raisimOgre raisim::assimp wrapperwalk2 nmpc_lib optimization_lib ${CASADI_LIBRARIES} ipopt)# osqp::osqpstatic)
target_include_directories(uprightwalk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(uprightwalk PRIVATE -std=c++17)

add_executable(uprightwalkR src/A1_up_walk_Raib.cpp)
target_link_libraries(uprightwalkR PUBLIC raisim::raisimOgre raisim::assimp wrapperwalk2 nmpcr_lib optimization_lib ${CASADI_LIBRARIES} ipopt)# osqp::osqpstatic)
target_include_directories(uprightwalkR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(uprightwalkR PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(uprightwalkR PRIVATE -std=c++17)

add_executable(uprightwalkExp src/A1_up_walk_Exp.cpp)
target_link_libraries(uprightwalkExp PUBLIC raisim::raisimOgre raisim::assimp wrapperwalk2 nmpc_lib optimization_lib ${CASADI_LIBRARIES} ipopt)# osqp::osqpstatic)
target_include_directories(uprightwalkExp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(uprightwalkExp PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(uprightwalkExp PRIVATE -std=c++17)

add_executable(makeUpright src/A1_makeupright.cpp)
target_link_libraries(makeUpright PUBLIC raisim::raisimOgre raisim::assimp wrapper2 nmpc_lib optimization_lib ${CASADI_LIBRARIES} ipopt)# osqp::osqpstatic)
target_include_directories(makeUpright PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(makeUpright PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(makeUpright PRIVATE -std=c++17)

add_executable(fullExpSim src/A1_Experiment_Script.cpp)
target_link_libraries(fullExpSim PUBLIC raisim::raisimOgre raisim::assimp wrapper2 nmpc_lib optimization_lib ${CASADI_LIBRARIES} ipopt)# osqp::osqpstatic)
target_include_directories(fullExpSim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(fullExpSim PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(fullExpSim PRIVATE -std=c++17)

add_executable(MT_Sim src/A1_Sim_MT.cpp)
target_link_libraries(MT_Sim PUBLIC raisim::raisimOgre raisim::assimp wrapper2 nmpc_lib optimization_lib ${CASADI_LIBRARIES} ${Boost_LIBRARIES} ipopt mscl OpenSSL::SSL OpenSSL::Crypto)# osqp::osqpstatic)
target_include_directories(MT_Sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS} ${MSCL_DIR}/MSCL/source)
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(MT_Sim PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
#target_include_directories(MT_Sim PRIVATE ${UT_INC})
target_compile_options(MT_Sim PRIVATE -std=c++17)