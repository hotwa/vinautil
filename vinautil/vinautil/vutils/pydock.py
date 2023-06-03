#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :pydock
@Description:       : pydock 调用相关类, 准备pydock配置文件
@Date               :2023/1/5 13:02:26
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional, NoReturn, List, Union, Dict, Tuple

'''
https://life.bsc.es/pid/pydock/doc/tutorial.html#setup-process
zdock配置参考：
对于蛋白小分子相互作用的配置参考：
在 [ligand] 部分中将 newmol 改为 B 是因为在 pyDock 输出文件中，配体的链名应与受体的链名不同。这是因为，在对齐过程中，配体可能会在受体上进行旋转和平移，并且在输出文件中，应为每个复合物配置提供单独的链名。因此，在这种情况下，将 newmol 改为 B 可确保在 pyDock 输出文件中，配体的链名与受体的链名不同。
如果不这样做，pyDock 可能无法正常运行。例如，如果在 [ligand] 部分中将 newmol 保留为 A，则在输出文件中受体和配体的链名将相同，这可能会导致 pyDock 无法识别受体和配体。
对接完成之后可以将小分子的链重新标记为A
[receptor]
pdb     = 1C5K.pdb
mol     = A
newmol  = A
restr   = A.His.246,A.Ala.249

[ligand]
pdb     = 1OAP.pdb
mol     = A
newmol  = B
restr   = B.Ala.88

对于蛋白蛋白相互作用的配置参考：
[receptor]
pdb		= 1C5K.pdb
mol		= A
newmol	= A
restr   = 

[ligand]
pdb		= 1OAP.pdb
mol		= B
newmol	= B
restr   = 
'''

@dataclass()
class pyDockBaseConfig():
    pdb: Path
    mol: str
    newmol: str
    restr: (str, type(None)) = field(default=None)

@dataclass
class pyDockConfigReceptor(pyDockBaseConfig):
    ...

@dataclass
class pyDockConfigLigand(pyDockBaseConfig):
    ...

@dataclass
class pyDockConfig():
    receptor: pyDockConfigReceptor
    ligand: pyDockConfigLigand

    def __str__(self):
        return f'''
pyDockConfig class
        '''

    def to_file(self, out_file: Path) -> NoReturn:
        with open(out_file.as_posix(), 'w', encoding='utf-8') as f:
            f.write(f'''
[receptor]
pdb = {self.receptor.pdb.as_posix()}
mol = {self.receptor.mol}
newmol  = {self.receptor.newmol}
restr   = {self.receptor.restr}

[ligand]
pdb = {self.ligand.pdb.as_posix()}
mol = {self.ligand.mol}
newmol  = {self.ligand.newmol}
restr   = {self.ligand.restr}
''')

    def to_dict(self):
        return asdict(self)


def testUnit():
    c = pyDockConfig(receptor=pyDockConfigReceptor(pdb=Path('312w.pdb'), mol='A', newmol='A', restr=None),
                 ligand=pyDockConfigLigand(pdb=Path('a12w.pdb'), mol='B', newmol='B', restr=None))
    c.to_file(Path('test.txt'))
    return c

