from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess

# Get GSL include and library paths using gsl-config
def get_gsl_flags(flag):
    return subprocess.check_output(["gsl-config", flag], text=True).strip().split()

ext_modules = [
    Pybind11Extension(
        'pcms.tree',
        [
            'src/tree.cpp',
            'src/dist.cpp',
            'src/pybind11.cpp'
        ],
        include_dirs=['include/'] + get_gsl_flags("--cflags"),
        library_dirs=get_gsl_flags("--libs") + get_gsl_flags("--prefix") + ["/usr/lib", "/usr/local/lib"],
        libraries=[
            "gsl", 
            "gslcblas"
        ],
        extra_compile_args=['-std=c++17'],
        extra_link_args=get_gsl_flags("--libs"),
    ),
]

setup(
    name='pcms',
    version="0.1.0",
    packages=find_packages(where='pcms'),
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
)