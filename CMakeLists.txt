cmake_minimum_required(VERSION 3.24)
project(DistributedComputing C)

find_package(MPI)

set(CMAKE_C_STANDARD 23)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(DistributedComputing main.c)
target_link_libraries(DistributedComputing ${MPI_C_LIBRARIES} dl m)
