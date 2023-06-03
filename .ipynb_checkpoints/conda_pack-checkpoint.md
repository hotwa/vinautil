## install vinautil

```shell
conda create -n my_env_test -c pylyzeng vinautil -y
```

## build anaconda package

```shell
conda create -n condabuild python=3 -y
conda activate condabuild
conda install conda-build anaconda-client -y
cd <source code dir>
anaconda login
conda build .
anaconda upload /path/to/conda-package.tar.bz2
```


