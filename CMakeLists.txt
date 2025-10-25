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

# --- Python and Pybind11 ---
find_package(
    Python 
    COMPONENTS Interpreter Development.Module
    REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# --- Use standard CMake build type flags ---
# Debug: no optimization, include debug symbols
# Release: full optimization
# RelWithDebInfo: optimization + debug info
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -fopenmp -funroll-loops -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -march=native -flto -fopenmp -funroll-loops -ffast-math")

message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

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

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "Enabling performance optimizations for ${CMAKE_BUILD_TYPE} build...")

    # Add compiler-specific optimization flags
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(_tree PRIVATE
            -O3                 # Full optimization
            -march=native       # Use all CPU features available on this machine
            -ffast-math         # Allow unsafe but faster math
            -funroll-loops      # Unroll loops for performance
            -fomit-frame-pointer
        )
    elseif(MSVC)
        target_compile_options(_tree PRIVATE
            /O2                 # Optimize for speed
            /fp:fast            # Fast math
            /openmp             # Enable OpenMP
        )
    endif()
endif()

# --- pcms._haar module ---
pybind11_add_module(_haar 
                    src/haar/pybind11.cpp
                    src/tree/tree.cpp 
                    src/tree/tree-dist.cpp
                    src/haar/haar-sparsify.cpp)
target_include_directories(_haar
                           PRIVATE 
                           ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(_haar 
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_haar
                      PRIVATE ${MKL_LIBS} ${SYSTEM_LIBS})
install(TARGETS _haar LIBRARY DESTINATION pcms)

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "Enabling performance optimizations for ${CMAKE_BUILD_TYPE} build...")

    # Try to enable OpenMP
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found — enabling parallelization")
        target_link_libraries(_haar PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP not found — continuing without parallelization")
    endif()

    # Add compiler-specific optimization flags
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(_haar PRIVATE
            -O3                 # Full optimization
            -march=native       # Use all CPU features available on this machine
            -ffast-math         # Allow unsafe but faster math
            -funroll-loops      # Unroll loops for performance
            -fomit-frame-pointer
        )
    elseif(MSVC)
        target_compile_options(_haar PRIVATE
            /O2                 # Optimize for speed
            /fp:fast            # Fast math
            /openmp             # Enable OpenMP
        )
    endif()
endif()