# install shell script for the project
$PYTHON -m pip install . --no-deps --ignore-installed -vv
git clone https://mirror.ghproxy.com/https://github.com/forlilab/Meeko
cd Meeko
$PYTHON -m pip install .  --no-deps --ignore-installed -vv