## install vinautil

```shell
conda create -n my_env_test -c pylyzeng -c conda-forge -c bioconda vinautil -y
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

```shell
micromamba create -n condabuild python=3 conda-build anaconda-client -c conda-forge -y
micromamba activate condabuild
cd <source code dir>
anaconda login
conda build .
anaconda upload /path/to/conda-package.tar.bz2
```


