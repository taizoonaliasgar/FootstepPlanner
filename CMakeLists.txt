cmake_minimum_required(VERSION 3.5)
project(footstep_planner)

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3")
set(USE_QPSWIFT TRUE) 
add_definitions(-DUSE_QPSWIFT=1)

find_package(raisim CONFIG HINTS "/home/taizoon/raisimEnv/raisimBuild/lib/cmake/raisim")
find_package(raisimOgre CONFIG HINTS "/home/taizoon/raisimEnv/raisimBuild/share/raisimOgre/cmake")
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



# ===== Parameter library ===== #
add_library(params_lib
    "${CMAKE_CURRENT_SOURCE_DIR}/params/Parameters.cpp")
target_include_directories(params_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/params
    ${CMAKE_CURRENT_SOURCE_DIR}/global_include)


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


add_executable(upright_NMPC src/upright_NMPC.cpp)
target_link_libraries(upright_NMPC PUBLIC raisim::raisimOgre raisim::assimp nmpc_lib ${CASADI_LIBRARIES} ipopt)
target_include_directories(upright_NMPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
target_include_directories(upright_NMPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(upright_NMPC PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(upright_NMPC PRIVATE -std=c++17)

add_executable(upright_NMPCR src/upright_NMPCR.cpp)
target_link_libraries(upright_NMPCR PUBLIC raisim::raisimOgre raisim::assimp nmpcr_lib ${CASADI_LIBRARIES} ipopt)
target_include_directories(upright_NMPCR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/raisimOgre_include ${EIGEN3_INCLUDE_DIRS})
target_include_directories(upright_NMPCR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(upright_NMPCR PRIVATE "-DEXAMPLE_ROBOT_RESOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/rsc/")
target_compile_options(upright_NMPCR PRIVATE -std=c++17)

