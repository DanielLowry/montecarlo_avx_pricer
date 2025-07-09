# setup.py
from skbuild import setup

setup(
    name="mc_pricer_py",
    version="0.1.0",
    packages=["mc_pricer_py"],         
    package_dir={"mc_pricer_py": "mc_pricer_py"},  
    cmake_source_dir="../analytics",
    cmake_args=["-DCMAKE_BUILD_TYPE=Release"],
    python_requires=">=3.7",
)