English | [中文](./README_cn.md)

## Description

All the codes in this repository were developed by Lingyu Zeng while he was pursuing a graduate degree in the School of Biological Engineering and Food at Hubei University of Technology. Some of the codes can be used to reproduce the results of Zeng's work during his master's degree, and the codes are for reference.

The code primarily focuses on the development of [AutoDock Vina](https://vina.scripps.edu/) and [SCARdock](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00334).

[SCARdock web site](https://scardock.com) for covalent inhibitor screening.

The package can also handle PDB files as well as small molecule files.

## Environment require

```shell
python >3.8
conda-forge::pandas >1.2.0
conda-forge::openbabel >3.0.0
conda-forge::rdkit
conda-forge::pymol-open-source
conda-forge::loguru
conda-forge::numpy
conda-forge::swig
conda-forge::boost-cpp
conda-forge::sphinx
conda-forge::sphinx_rtd_theme
conda-forge::vina
conda-forge::ipython
conda-forge::peewee
conda-forge::scipy
conda-forge::prody  # meeko dependency (optionally, for covalent docking)
```

[more package](./environment.yml)

## Quick Start

### Install

```shell
conda create -n vina_env python=3 --yes
conda activate vina_env
conda install -c pylyzeng vinautil
```

### Build method

```shell
python .\setup.py sdist # build vinautil.tar.gz
conda build -c conda-forge . # build conda package
```

## Acknowledgement

[Sen Liu](https://sgsp.hbut.edu.cn/info/1085/1794.htm)

[Qi Song](https://sgsp.hbut.edu.cn/info/1087/1813.htm)

[Lab site](http://www.liugroup.site)

## License
[MIT](./LICENSE)