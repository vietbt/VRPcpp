import os
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

def load_cpp_files(folder):
    files = [os.path.join(folder, file) for file in os.listdir(folder) if file.endswith(".cpp")]
    return files

__version__ = "1.0.3"

ext_modules = [
    Pybind11Extension(
        "evrp_cpp",
        load_cpp_files("VRPcpp/vns"),
        define_macros=[('VERSION_INFO', __version__)],
    ),
    Pybind11Extension(
        "cvrp_cpp",
        load_cpp_files("VRPcpp/hgs"),
        define_macros=[('VERSION_INFO', __version__)],
    ),
]

setup(
    name="vrp_cpp",
    version=__version__,
    author="VietBT",
    author_email="tvbui@smu.edu.sg",
    description="VRP C++ Wrapper",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=True,
    python_requires=">=3.6",
)
