cmake_minimum_required(VERSION 3.12)
project(openvdb-remesher)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenVDB REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(EIGEN3 REQUIRED NO_MODULE)

set(VCGLIB_DIR "../.." CACHE INTERNAL "Path to vcglib")
set(EIGEN_DIR  = "${VCGLIB_DIR}/eigenlib" CACHE INTERNAL "")


add_executable(${PROJECT_NAME} OpenVDBRemeshing.cpp
    ${VCGLIB_DIR}/wrap/ply/plylib.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${VCGLIB_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${TBB_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${OpenVDB_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN_DIR})


target_link_libraries(${PROJECT_NAME} PRIVATE
  ${OpenVDB_LIBRARIES}
  tbb
  ${Boost_LIBRARIES}
)