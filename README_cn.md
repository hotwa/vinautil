English(./README.md) | [中文]

## 简介
本仓库所有的代码，为曾令宇在湖北工业大学生物工程与食品学院攻读研究生学位时开发。部分代码可以用于复现曾令宇硕士期间的工作成果，代码供参考。

这些代码主要是为 [autodock vina](https://vina.scripps.edu/) 和 [SCARdock](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00334) 开发。

[SCARdock 网站](https://scardock.com)筛选共价抑制剂。

同时也可以处理PDB文件以及小分子文件

## 环境

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
- conda-forge::prody # meeko dependecy (optionally, for covalent docking)
```

## 使用

### 安装

```shell
conda create -n vinautil_env -c pylyzeng vinautil --yes
conda activate vinautil_env
```

### 测试

```shell
scardocktest
```

### 对接

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

### 代码安装

推荐使用conda安装

### [构建](./conda_pack.md)

## 致谢

[Sen Liu](https://sgsp.hbut.edu.cn/info/1085/1794.htm)

[Qi Song](https://sgsp.hbut.edu.cn/info/1087/1813.htm)

[Lab site](http://www.liugroup.site)

## 许可协议
[MIT](./LICENSE)