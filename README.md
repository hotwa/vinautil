English | [中文](./README_cn.md)

## Description

All the codes in this repository were developed by Lingyu Zeng while he was pursuing a graduate degree in the School of Biological Engineering and Food at Hubei University of Technology. Some of the codes can be used to reproduce the results of Zeng's work during his master's degree, and the codes are for reference.

The code is mainly developed for [autodock vina](https://vina.scripps.edu/) and [SCARdock](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00334)

The package can also handle PDB files as well as small molecule files.

## Environment require

```shell
python >3.9
openbabel >3.0.0
pandas >1.2.0
rdkit =2022.03.5
pymol-open-source =2.5.0
loguru
```

[more package](./environment.yml)

## Quick Start

### Install

```shell
conda create -n vina_env python=3.9 --yes
conda activate vina_env
conda install -c pylyzeng vinautil
```

### Install from source

```shell
git clone https://github.com/hotwa/vinautil.git
cd vinautil
python setup.py install
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