from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess
import os

DEBUG = os.getenv("DEBUG_BUILD", "0") == "1"

extra_compile_args = ['-std=c++17']
if DEBUG:
    print("DEBUG_BUILD set")
    extra_compile_args += ['-g', '-O0']
else:
    extra_compile_args += ['-O3']

def get_gsl_flags(flag):
    return subprocess.check_output(['gsl-config', flag], text=True).strip().split()

ext_modules = [
    Pybind11Extension(
        'pcms._tree',
        sources=[
            'src/tree/tree.cpp',
            'src/tree/tree-dist.cpp',
            'src/tree/pybind11.cpp'
        ],
        include_dirs=['include/'],
        library_dirs=get_gsl_flags('--libdir') + ['/usr/lib', '/usr/local/lib'],
        libraries=['gsl', 'gslcblas', 'm'],
        extra_compile_args=extra_compile_args,
        extra_link_args=['-g'],
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
        extra_compile_args=extra_compile_args,
        extra_link_args=['-g'],
        language='c++',
        define_macros=[('DEBUG', '1')] if DEBUG else []
    )
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
