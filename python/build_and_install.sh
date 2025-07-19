
# activate the venv
source ./venv/bin/activate

# build the analytics, and package into a wheel using skbuild
python3 setup.py bdist_wheel

# install to the venv - overwrite any previous package
# currently not very portable..
pip install ./dist/mc_pricer_py-0.1.0-cp312-cp312-linux_x86_64.whl --force-reinstall