# install shell script for the project
$PYTHON -m pip install . --no-deps --ignore-installed -vv
git clone https://github.com/forlilab/Meeko
cd Meeko
pip install .  --no-deps --ignore-installed -vv
#$PYTHON setup.py install