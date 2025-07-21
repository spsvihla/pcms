from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess
import os

DEBUG = os.getenv("DEBUG_BUILD", "0") == "1"

# Common compiler and linker flags
common_compile_args = ['-std=c++17']
common_link_args = []

if DEBUG:
    print("DEBUG_BUILD set")
    common_compile_args += ['-g', '-O0', '-march=native']
    common_link_args += ['-g']
else:
    common_compile_args += ['-O3', '-march=native', '-flto', '-fopenmp', '-ffast-math']
    common_link_args += ['-flto', '-fopenmp']

def get_gsl_lib_dirs():
    try:
        return subprocess.check_output(['gsl-config', '--libdir'], text=True).strip().split()
    except Exception:
        return []

gsl_lib_dirs = get_gsl_lib_dirs() + ['/usr/lib', '/usr/local/lib']

ext_modules = [
    Pybind11Extension(
        'pcms._tree',
        sources=[
            'src/tree/tree.cpp',
            'src/tree/tree-dist.cpp',
            'src/tree/pybind11.cpp'
        ],
        include_dirs=['include/'],
        library_dirs=gsl_lib_dirs,
        libraries=['gsl', 'gslcblas', 'm'],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        language='c++',
        define_macros=[('DEBUG', '1')] if DEBUG else []
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
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        language='c++',
        define_macros=[('DEBUG', '1')] if DEBUG else []
    ),
]

setup(
    name='pcms',
    version='0.1.0',
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
