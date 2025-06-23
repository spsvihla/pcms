from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess, os


DEBUG = os.getenv("DEBUG_BUILD", "0") == "1"

# Set compile args depending on debug mode 
extra_compile_args=['-std=c++17']
if DEBUG:
    print("DEBUG_BUILD set")
    extra_compile_args += ["-g","-O0"]
else:
    extra_compile_args += ["-O3"]

# Get GSL include and library paths using gsl-config
def get_gsl_flags(flag):
    return subprocess.check_output(["gsl-config", flag], text=True).strip().split()

ext_modules = [
    Pybind11Extension(
        'pcms._tree',
        [
            'src/tree/tree.cpp',
            'src/tree/tree-dist.cpp',
            'src/tree/pybind11.cpp'
        ],
        include_dirs=['include/'],
        library_dirs=get_gsl_flags("--libs") + get_gsl_flags("--prefix") + ["/usr/lib", "/usr/local/lib"],
        libraries=[
            "gsl", 
            "gslcblas"
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=get_gsl_flags("--libs"),
    ),
    Pybind11Extension(
        'pcms._haar',
        sources=[
            'src/tree/tree.cpp',
            'src/haar/haar-dist.cpp',
            'src/haar/haar-sparsify.cpp',
            'src/haar/pybind11.cpp'
        ],
        include_dirs=['include/'],
        extra_compile_args=extra_compile_args       
    )
]

setup(
    name='pcms',
    version="0.1.0",
    packages=find_packages(),
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    extras_require={
        'notebooks': [
            'numpy',
            'scipy',
            'pandas',
            'matplotlib',
            'jupyter'
        ],
    },
)