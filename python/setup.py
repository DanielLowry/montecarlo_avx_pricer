# setup.py
from skbuild import setup
import os

# Check for an environment variable or command-line flag
enable_logging = os.environ.get("ENABLE_LOGGING", "OFF")

setup(
    name="mc_pricer_py",
    version="0.1.0",
    packages=["mc_pricer_py"],         
    cmake_source_dir="../analytics",
    cmake_args=[
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DENABLE_LOGGING={enable_logging}"
    ],
    python_requires=">=3.7",
)