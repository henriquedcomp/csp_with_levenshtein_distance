cmake_minimum_required(VERSION 3.10.0)
project(csp_with_levenshtein_distance VERSION 0.1.3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#OR-Tools config
find_package(ortools CONFIG REQUIRED)

#Main executable
add_executable(
    csp_with_levenshtein_distance
    src/main.cpp

    src/hayashida_and_koyano_2017_mpsolver.cpp
    src/hayashida_and_koyano_2017_cpsolver.cpp
    src/arbib_et_al_binary_path_mpsolver.cpp
    src/arbib_et_al_binary_path_cpsolver.cpp
    src/arbib_et_al_path_mpsolver.cpp
    src/improved_hayashida_and_koyano_2017_mpsolver.cpp
    src/improved_hayashida_and_koyano_2017_cpsolver.cpp
    src/improved_arbib_et_al_path_cpsolver.cpp
    src/improved_arbib_et_al_path_mpsolver.cpp
    src/improved_arbib_et_al_binary_path_mpsolver.cpp
    src/improved_arbib_et_al_binary_path_cpsolver.cpp
    src/arbib_et_al_path_cpsolver.cpp
)
target_include_directories(csp_with_levenshtein_distance PRIVATE include)
target_link_libraries(csp_with_levenshtein_distance PRIVATE ortools::ortools)

#Instance generator
add_executable(
    instance_generator
    instance_generator.cpp
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_compile_options(csp_with_levenshtein_distance PRIVATE 
    $<$<CXX_COMPILER_ID:GNU>:
        -O3 
        -march=native 
        -flto=auto
    >
)

target_link_options(csp_with_levenshtein_distance PRIVATE 
    $<$<CXX_COMPILER_ID:GNU>:
        -flto
    >
)
