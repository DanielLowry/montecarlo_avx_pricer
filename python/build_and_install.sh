
# Usage: ./build_and_install.sh [ON|OFF]
LOGGING=${1:-OFF}  # Default to OFF if not provided


# activate the venv
source ./venv/bin/activate

# build the analytics, and package into a wheel using skbuild
ENABLE_LOGGING=$LOGGING python3 setup.py bdist_wheel

# install to the venv - overwrite any previous package
# currently not very portable..
pip install ./dist/mc_pricer_py-0.1.0-cp312-cp312-linux_x86_64.whl --force-reinstall