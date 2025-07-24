import sys
import subprocess
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


def parse_args():
    """
    Parse --build-type={debug|profile|release} (default: release).
    Removes it from sys.argv for setuptools.
    """
    build_type = "release"
    filtered_args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg.startswith("--build-type="):
            build_type = arg.split("=", 1)[1].lower()
            if build_type not in {"debug", "profile", "release"}:
                raise ValueError(f"Unknown build type: {build_type}")
            i += 1
        else:
            filtered_args.append(arg)
            i += 1

    sys.argv = [sys.argv[0]] + filtered_args
    return build_type


def get_lib_flags(build_type: str) -> tuple[list[str], list[str], list[str]]:
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
    mkl_includes = ["/opt/intel/oneapi/mkl/latest/include"]
    mkl_libdirs = ["/opt/intel/oneapi/mkl/latest/lib/intel64"]
    mkl_libs = ["mkl_intel_lp64", "mkl_core"]
    if build_type in ["release", "profile"]:
        mkl_libs += ["mkl_intel_thread", "iomp5"]       # multi-threaded
    else:
        mkl_libs += ["mkl_sequential"]                  # single-threaded

    sys_libs = ["m"]

    include_dirs = gsl_includes + mkl_includes
    library_dirs = ["/usr/lib", "/usr/local/lib"] + gsl_libdirs + mkl_libdirs
    libraries = gsl_libs + mkl_libs + sys_libs

    return include_dirs, library_dirs, libraries


def get_compile_and_link_args(build_type: str) -> tuple[list[str], list[str]]:
    compile_args = ["-std=c++20"]
    link_args = []

    if build_type == "debug":
        print("Build type: DEBUG")
        compile_args.extend(["-g", "-O0"])
        link_args.append("-g")
    elif build_type == "profile":
        print("Build type: PROFILE")
        compile_args.extend(["-g", "-O3", "-march=native", "-flto", "-fopenmp", "-ffast-math"])
        link_args.extend(["-g", "-flto", "-fopenmp"])
    elif build_type == "release":
        print("Build type: RELEASE")
        compile_args.extend(["-O3", "-march=native", "-flto", "-fopenmp", "-ffast-math"])
        link_args.extend(["-flto", "-fopenmp"])
    else:
        raise ValueError(f"Unknown build type: {build_type}")

    return compile_args, link_args


def main():
    build_type = parse_args()

    compile_args, link_args = get_compile_and_link_args(build_type)

    include_dirs, library_dirs, libraries = get_lib_flags(build_type)

    debug_macro = [("DEBUG", "1")] if build_type in {"debug", "profile"} else []

    ext_modules = [
        Pybind11Extension(
            "pcms._tree",
            sources=[
                "src/tree/tree.cpp",
                "src/tree/tree-dist.cpp",
                "src/tree/pybind11.cpp"
            ],
            include_dirs=["include/"] + include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            language="c++",
            define_macros=debug_macro,
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
            define_macros=debug_macro,
        ),
    ]

    setup(
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext},
    )


if __name__ == "__main__":
    main()
