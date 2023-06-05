English | [中文](./README_cn.md)

## Introduction

All the code in this repository was developed by Zeng Lingyu during his pursuit of a graduate degree at the College of Biological Engineering and Food Science, Hubei University of Technology. Some of the code can be used to reproduce Zeng Lingyu's research work during his master's degree and is provided for reference.

These codes are mainly developed for [AutoDock Vina](https://vina.scripps.edu/) and [SCARdock](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00334).

[SCARdock website](https://scardock.com/) for screening covalent inhibitors.

It can also handle PDB files and small molecule files.

## Environment require

```yaml
- python =3.10
- bioconda::mgltools
- conda-forge::spyrmsd
- conda-forge::pandas
- conda-forge::openbabel
- conda-forge::rdkit
- conda-forge::pymol-open-source 2.5.0 py310h292d129_6
- conda-forge::loguru
- conda-forge::swig
- conda-forge::boost-cpp 
- conda-forge::sphinx
- conda-forge::sphinx_rtd_theme
- conda-forge::vina 1.2.3 py310he924329_2
- conda-forge::ipython
- conda-forge::biopython
- conda-forge::prody # meeko dependency (optionally, for covalent docking)
```

## Usage

### Installation

```shell
conda create -n vinautil_env -c pylyzeng vinautil --yes
conda activate vinautil_env
```

### Test

```shell
# before test need execute:
scardock -h
scardocktest
```

### Docking

```shell
scardock --help
usage: scardock [-h] [-r [recepotr file]] [-l [ligand file]] [-s [residue covalent site]] [-c [covalent chain ID]]
                [-log [output log directory]]

SCARdock Docking

options:
  -h, --help            show this help message and exit
  -r [recepotr file], --receptor [recepotr file]
                        recepotr file, support pdb
  -l [ligand file], --ligand [ligand file]
                        ligand file, molecule file (MOL2, SDF,...)(use meeko prepare)
  -s [residue covalent site], --site [residue covalent site]
                        residue covalent site
  -c [covalent chain ID], --chain [covalent chain ID]
                        covalent chain ID
  -log [output log directory], --log_dir [output log directory]
                        Relative Path
```

## Code Installation

It is recommended to use conda for installation.

## Acknowledgement

[Sen Liu](https://sgsp.hbut.edu.cn/info/1085/1794.htm)

[Qi Song](https://sgsp.hbut.edu.cn/info/1087/1813.htm)

[Lab site](http://www.liugroup.site)

## License
[MIT](./LICENSE)