English(./README.md) | [中文]

## 简介
本仓库所有的代码，为曾令宇在湖北工业大学生物工程与食品学院攻读研究生学位时开发。部分代码可以用于复现曾令宇硕士期间的工作成果，代码供参考。

这些代码主要是为 [autodock vina](https://vina.scripps.edu/) 和 [SCARdock](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00334) 开发

同时也可以处理PDB文件以及小分子文件

## 环境

```shell
python >3.9
openbabel >3.0.0
pandas >1.2.0
rdkit =2022.03.5
pymol-open-source =2.5.0
loguru
```

[更多参考](./environment.yml)

## 使用

### 安装

```shell
conda create -n vina_env python=3.9 --yes
conda activate vina_env
conda install -c pylyzeng vinautil
```

### 代码安装

```shell
git clone https://github.com/hotwa/vinautil.git
cd vinautil
python setup.py install
```

### 构建

```shell
python .\setup.py sdist # build vinautil.tar.gz
conda build -c conda-forge . # build conda package
```

## 致谢

[Sen Liu](https://sgsp.hbut.edu.cn/info/1085/1794.htm)

[Qi Song](https://sgsp.hbut.edu.cn/info/1087/1813.htm)

[Lab site](http://www.liugroup.site)

## 许可协议
[MIT](./LICENSE)