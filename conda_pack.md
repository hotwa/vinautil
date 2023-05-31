before use this package：

## install autodock vina 1.2.3

```shell
$ conda create -n vina python=3
$ conda activate vina
$ conda config --env --add channels conda-forge
$ conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
$ pip install vina
```

test code

```python
from Vina import vina
```

