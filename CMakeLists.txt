cmake_minimum_required(VERSION 3.15...3.26)
project(
    "${SKBUILD_PROJECT_NAME}"
    LANGUAGES CXX
    VERSION "${SKBUILD_PROJECT_VERSION}")

# --- Parse build type argument
if(NOT DEFINED BUILD_TYPE)
    set(BUILD_TYPE "Release")
endif()
string(TOUPPER "${BUILD_TYPE}" BUILD_TYPE_UPPER)

message(STATUS "Build type: ${BUILD_TYPE_UPPER}")

# --- Set C++ standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Pythong and Pybind11 ---
find_package(
    Python 
    COMPONENTS Interpreter Development.Module
    REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# --- Compiler flags per build type ---
if(BUILD_TYPE_UPPER STREQUAL "DEBUG")
    add_compile_options(-g -O0)
    add_link_options(-g)
elseif(BUILD_TYPE_UPPER STREQUAL "PROFILE")
    add_compile_options(-g -O3 -march=native -flto -fopenmp -funroll-loops -ffast-math)
    add_link_options(-g -flto -fopenmp)
elseif(BUILD_TYPE_UPPER STREQUAL "RELEASE")
    add_compile_options(-O3 -march=native -flto -fopenmp -funroll-loops -ffast-math)
    add_link_options(-flto -fopenmp)
else()
    message(FATAL_ERROR "Unknown BUILD_TYPE: ${BUILD_TYPE}")
endif()

# --- MKL & system libs ---
set(MKL_INCLUDE_DIRS "/opt/intel/oneapi/mkl/latest/include")
set(MKL_LIBRARY_DIRS "/opt/intel/oneapi/mkl/latest/lib/intel64")
set(MKL_LIBS "mkl_intel_lp64" "mkl_core")
if(BUILD_TYPE_UPPER STREQUAL "DEBUG")
    list(APPEND MKL_LIBS "mkl_sequential")              # single-threaded
else()
    list(APPEND MKL_LIBS "mkl_intel_thread" "iomp5")    # multi-threaded
endif()
set(SYSTEM_LIBS "m")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${MKL_INCLUDE_DIRS})
link_directories(${MKL_LIBRARY_DIRS} /usr/lib /usr/local/lib)

# --- pcms._tree module ---
pybind11_add_module(_tree 
                    src/tree/pybind11.cpp
                    src/tree/tree.cpp 
                    src/tree/tree-dist.cpp)
target_include_directories(_tree
                           PRIVATE 
                           ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(_tree 
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_tree
                      PRIVATE ${MKL_LIBS} ${SYSTEM_LIBS})
install(TARGETS _tree LIBRARY DESTINATION pcms)

# --- pcms._haar module ---
pybind11_add_module(_haar 
                    src/haar/pybind11.cpp
                    src/tree/tree.cpp 
                    src/tree/tree-dist.cpp
                    src/haar/haar-sparsify.cpp
                    src/haar/haar-dist.cpp)
target_include_directories(_haar
                           PRIVATE 
                           ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(_haar 
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_haar
                      PRIVATE ${MKL_LIBS} ${SYSTEM_LIBS})
install(TARGETS _haar LIBRARY DESTINATION pcms)
