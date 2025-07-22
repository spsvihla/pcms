import sys
import subprocess
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


def parse_args():
    """Parse custom command-line args and remove them from sys.argv."""
    debug = False
    optimize = False
    opt_level = None
    filtered_args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--debug":
            debug = True
            i += 1
        elif arg == "--optimize":
            optimize = True
            i += 1
        elif arg.startswith("--opt="):
            opt_level = arg.split("=", 1)[1]
            i += 1
        elif arg == "--opt":
            if i + 1 >= len(sys.argv):
                raise RuntimeError("--opt requires a value (e.g., --opt -O2)")
            opt_level = sys.argv[i + 1]
            i += 2
        else:
            filtered_args.append(arg)
            i += 1

    sys.argv = [sys.argv[0]] + filtered_args
    return debug, optimize, opt_level


def get_lib_flags() -> tuple[list[str], list[str], list[str]]:
    """
    Get combined include_dirs, library_dirs, and libraries for GSL + MKL.
    Adjust paths to match your system!
    """
    # ---- GSL ----
    try:
        gsl_libdir = subprocess.check_output(["gsl-config", "--libdir"], text=True).strip()
        gsl_includedir = subprocess.check_output(["gsl-config", "--cflags"], text=True).strip()
        gsl_includes = [flag[2:] for flag in gsl_includedir.split() if flag.startswith("-I")]
        gsl_libdirs = gsl_libdir.split()
        gsl_libs = ["gsl", "gslcblas"]
    except Exception as e:
        raise RuntimeError("Could not find GSL using gsl-config") from e

    # ---- MKL ----
    # Example for Intel oneAPI MKL
    mkl_includes = ["/opt/intel/oneapi/mkl/latest/include"]
    mkl_libdirs = ["/opt/intel/oneapi/mkl/latest/lib/intel64"]

    # Common static link MKL
    mkl_libs = ["mkl_intel_lp64", "mkl_core", "mkl_sequential"]

    # Always need system math lib too
    sys_libs = ["m"]

    include_dirs = gsl_includes + mkl_includes
    library_dirs = ["/usr/lib", "/usr/local/lib"] + gsl_libdirs + mkl_libdirs
    libraries = gsl_libs + mkl_libs + sys_libs

    return include_dirs, library_dirs, libraries


def get_compile_and_link_args(debug: bool, optimize: bool, opt_level: str | None) -> tuple[list[str], list[str]]:
    compile_args = ["-std=c++17"]
    link_args = []

    if debug:
        print("DEBUG_BUILD set: Using debug compile flags.")
        compile_args.append("-g")
        link_args.append("-g")

    if optimize:
        compile_args.extend(["-march=native", "-flto", "-fopenmp", "-ffast-math"])
        link_args.extend(["-flto", "-fopenmp"])

    if opt_level:
        compile_args.append(opt_level)
    elif optimize:
        compile_args.append("-O3")
    else:
        compile_args.append("-O0")

    return compile_args, link_args


def main():
    debug, optimize, opt_level = parse_args()

    compile_args, link_args = get_compile_and_link_args(debug, optimize, opt_level)

    include_dirs, library_dirs, libraries = get_lib_flags()

    ext_modules = [
        Pybind11Extension(
            "pcms._tree",
            sources=[
                "src/tree/tree.cpp",
                "src/tree/tree-dist.cpp",
                "src/tree/pybind11.cpp",
            ],
            include_dirs=["include/"] + include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            language="c++",
            define_macros=[("DEBUG", "1")] if debug else [],
        ),
        Pybind11Extension(
            "pcms._haar",
            sources=[
                "src/tree/tree.cpp",
                "src/haar/haar-dist.cpp",
                "src/haar/haar-sparsify.cpp",
                "src/haar/pybind11.cpp",
            ],
            include_dirs=["include/"] + include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            language="c++",
            define_macros=[("DEBUG", "1")] if debug else [],
        ),
    ]

    setup(
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext},
    )


if __name__ == "__main__":
    main()
