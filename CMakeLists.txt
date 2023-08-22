cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(dune2d LANGUAGES CXX C)
set(CMAKE_C_FLAGS $ENV{CFLAGS})
set(CMAKE_CXX_STANDARD 17)
find_package(ROOT REQUIRED)

message("status of ROOT_FOUND ${ROOT_FOUND}")

set(DK2NU_INCLUDE_DIR $ENV{DK2NU_INC})
set(DK2NU_LIB_DIR $ENV{DK2NU_LIB})

include_directories(SYSTEM
            ${ROOT_INCLUDE_DIRS}
            ${DK2NU_INCLUDE_DIR}
)

set(libraries ${libraries}
              ${ROOT_LIBRARIES}
              ${DK2NU_LIB_DIR}/libdk2nuTree.so
             )

add_executable(dune2Dflux
               extract_duneflux.cpp
              dune2Dflux.cpp
)
target_link_libraries(dune2Dflux ${libraries})