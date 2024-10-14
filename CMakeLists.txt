cmake_minimum_required(VERSION 3.12)
project(footstep_planner)

# if(UNIX AND NOT APPLE)
#     set(CMAKE_C_COMPILER "/usr/bin/gcc")
#     set(CMAKE_CXX_COMPILER "/usr/bin/g++")
# elseif(APPLE)
#     set(CMAKE_C_COMPILER "clang")
#     set(CMAKE_CXX_COMPILER "clang++")
# else()
#     message(WARNING "Only support LINUX and APPLE for now.")
# endif()

# for color output only
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

SET(CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,${CMAKE_SOURCE_DIR}/build")

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-O3")
set(USE_QPSWIFT TRUE) 
add_definitions(-DUSE_QPSWIFT=1)

# set the g++ compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")


if(CMAKE_CXX_COMPILER_LOADED)
    message(STATUS "Compiler Path: ${CMAKE_CXX_COMPILER}")
    message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
endif()

# set(LL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../A1_LL_Exp-Arch_Change)
set(LL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/LL_integrator)
set(ROBOT_DIR "A1_New_wFoot")

#find_package(raisim CONFIG HINTS "/home/taizoon/raisimEnv/raisimBuild/lib/cmake/raisim")
#find_package(raisimOgre CONFIG HINTS "/home/taizoon/raisimEnv/raisimBuild/share/raisimOgre/cmake")

find_package(raisim CONFIG HINTS "/home/trec/WorkRaj/raisim_legged/raisim_build/lib/cmake/raisim")
find_package(OGRE CONFIG REQUIRED PATHS "/home/trec/WorkRaj/raisim_legged/raisim_build/lib/OGRE/cmake")
find_package(raisimOgre CONFIG HINTS  "/home/trec/WorkRaj/raisim_legged/raisim_build/share/raisimOgre/cmake")  

include_directories(PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include)
find_package(Eigen3 REQUIRED HINTS ${Eigen3_HINT})

#find_package(osqp REQUIRED)

find_package(casadi REQUIRED)
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

include_directories(PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

set(dyn_srcs 
    "${LL_DIR}/dynamics/A1/A1_Dynamics_upr.cpp"
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
    dyn_lib
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
add_library(wrapperwalk
    "${CMAKE_CURRENT_SOURCE_DIR}/src/LocoWrapperwalk.cpp")
target_include_directories(wrapperwalk PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
 
 target_link_libraries(wrapperwalk PUBLIC
    ctrl_lib
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

# add_library(nmpcr_lib
#     "src/SRBNMPCR.cpp")
# target_compile_options(nmpcr_lib PRIVATE -std=c++17)

# target_include_directories(nmpcr_lib PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
#     ${CMAKE_CURRENT_SOURCE_DIR}/util_include
#     ${CMAKE_CURRENT_SOURCE_DIR}/global_include)
# target_link_libraries(nmpcr_lib PUBLIC 
#     params_lib
#     -L/usr/local/lib)



add_executable(upright_NMPC src/upright_NMPC.cpp)
target_link_libraries(upright_NMPC PUBLIC raisim::raisimOgre raisim::assimp nmpc_lib ${CASADI_LIBRARIES} ipopt stdc++fs)
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
target_link_libraries(uprightwalk PUBLIC raisim::raisimOgre raisim::assimp wrapperwalk nmpc_lib optimization_lib ${CASADI_LIBRARIES} ipopt stdc++fs)# osqp::osqpstatic)
target_include_directories(uprightwalk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
#target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${LL_DIR}/rsc/")
target_compile_definitions(uprightwalk PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(uprightwalk PRIVATE -std=c++17)
