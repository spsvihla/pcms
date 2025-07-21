import sys
import subprocess
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


def parse_args():
    """Parse custom command-line args and remove them from sys.argv."""
    debug = False
    filtered_args = []
    for arg in sys.argv[1:]:
        if arg == "--debug":
            debug = True
        else:
            filtered_args.append(arg)
    # Replace sys.argv with filtered args to not confuse setuptools
    sys.argv = [sys.argv[0]] + filtered_args
    return debug


def get_gsl_lib_dirs() -> list[str]:
    try:
        output = subprocess.check_output(["gsl-config", "--libdir"], text=True).strip()
        return output.split()
    except Exception as e:
        raise RuntimeError("GSL not found. Please install GSL.") from e


def get_eigen_include() -> list[str]:
    try:
        output = subprocess.check_output(["pkg-config", "--cflags", "eigen3"], text=True).strip()
        flags = output.split()
        includes = [flag[2:] for flag in flags if flag.startswith("-I")]
        if not includes:
            raise RuntimeError("No Eigen include path found in pkg-config output.")
        return includes
    except Exception as e:
        raise RuntimeError("Eigen not found. Please install Eigen development files.") from e


def get_compile_and_link_args(debug: bool) -> tuple[list[str], list[str]]:
    if debug:
        print("DEBUG_BUILD set: Using debug compile flags.")
        compile_args = ["-std=c++17", "-g", "-O0", "-march=native"]
        link_args = ["-g"]
    else:
        compile_args = ["-std=c++17", "-O3", "-march=native", "-flto", "-fopenmp", "-ffast-math"]
        link_args = ["-flto", "-fopenmp"]
    return compile_args, link_args


def main():
    debug = parse_args()

    compile_args, link_args = get_compile_and_link_args(debug)

    gsl_lib_dirs = get_gsl_lib_dirs()
    eigen_include_dirs = get_eigen_include()

    ext_modules = [
        Pybind11Extension(
            "pcms._tree",
            sources=[
                "src/tree/tree.cpp",
                "src/tree/tree-dist.cpp",
                "src/tree/pybind11.cpp",
            ],
            include_dirs=["include/"],
            library_dirs=["/usr/lib", "/usr/local/lib"] + gsl_lib_dirs,
            libraries=["gsl", "gslcblas", "m"],
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
            include_dirs=["include/"] + eigen_include_dirs,
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
