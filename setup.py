from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        'pcms.tree',
        [
            'src/tree.cpp',
            'src/dist.cpp',
            'src/pybind11.cpp'
        ],
        include_dirs=['include/'],
        extra_compile_args=['-std=c++17'],
    ),
]

setup(
    name='pcms',
    version="0.1.0",
    packages=find_packages(where='pcms'),
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
)