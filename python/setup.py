# setup.py
from skbuild import setup

setup(
    name="mc_pricer",
    version="0.1.0",
    packages=["py_pricer"],
    include_package_data=True,
    package_data={"py_pricer" : ["py_pricer.so"]},
    cmake_source_dir="../analytics",
    cmake_args=["-DCMAKE_BUILD_TYPE=Release"],
    python_requires=">=3.7",
)