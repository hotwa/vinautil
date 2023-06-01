#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :vinaconfig.py
@Description:       :
@Date     :2022/11/11 15:17:33
@Author      :hotwa
@version      :1.0
'''
from pathlib import Path, PurePath
from dataclasses import dataclass
from typing import Iterable
from collections.abc import Iterable as _Iterator


from vinautil.utils.typecheck import typeassert


@typeassert(x=float, y=float, z=float)
@dataclass
class xyz_point:
    __slots__ = ['x', 'y', 'z', '__dict__']
    x: float
    y: float
    z: float

    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self):
        return 'xyz point: {},{},{}'.format(self.x, self.y, self.z)

    @classmethod
    def from_iterable(cls, iterable_object: Iterable):
        if isinstance(iterable_object, _Iterator):
            return cls(*iterable_object)

@typeassert(receptor=PurePath,
            ligand=PurePath,
            box_size=xyz_point,
            center=xyz_point,
            outpdbqtfile=PurePath,
            exhaustiveness=int,
            num_modes=int,
            energy_range=int)
@dataclass
class vinaconfig:
    receptor: PurePath
    ligand: PurePath
    box_size: xyz_point
    center: xyz_point
    outpdbqtfile: PurePath
    exhaustiveness: int
    num_modes: int
    energy_range: int

    __slots__ = ['receptor', 'ligand', 'center',
                 'box_size', 'outpdbqtfile', '__dict__',
                 'exhaustiveness', 'num_modes', 'energy_range']

    def __init__(self,
                 receptor, ligand,
                 center,
                 box_size,
                 outpdbqtfile,
                 exhaustiveness = 32,
                 num_modes = 20,
                 energy_range = 5,
                 ):
        self.receptor = receptor
        self.ligand = ligand
        self.center = center if isinstance(center, xyz_point) else xyz_point.from_iterable(center)
        self.box_size = box_size if isinstance(box_size, xyz_point) else xyz_point.from_iterable(box_size)
        self.outpdbqtfile = outpdbqtfile
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.energy_range = energy_range

    def to_dict(self )->dict:
        return {'receptor': self.receptor,
                'ligand': self.ligand,
                'center': self.center,
                'box_size': self.box_size,
                'outpdbqtfile': self.outpdbqtfile,
                'exhaustiveness': self.exhaustiveness,
                'num_modes': self.num_modes,
                'energy_range': self.energy_range}

    def to_txt(self, file :Path,
               cx = None,
               cy = None,
               cz = None,
               sx = None,
               sy = None,
               sz = None,
               exhaustiveness = None,
               num_modes = None,
               energy_range = None,
               absolute_path=False) -> None:
        """to_txt save config to file supported autodock vina 1.1.2 and 1.2.3 binary file

        use: vina --config config.txt

        Arguments:
            file {file path} -- out put file
        """
        if absolute_path:
            receptor_path = self.receptor.absolute().as_posix()
            ligand_path = self.ligand.absolute().as_posix()
            outpdbqtfile_path = self.outpdbqtfile.absolute().as_posix()
        else:
            receptor_path = self.receptor.as_posix()
            ligand_path = self.ligand.as_posix()
            outpdbqtfile_path = self.outpdbqtfile.as_posix()
        exhaustiveness = exhaustiveness if exhaustiveness else self.exhaustiveness
        num_modes = num_modes if num_modes else self.num_modes
        energy_range = energy_range if energy_range else self.energy_range
        cx = cx if cx else self.center.x
        cy = cy if cy else self.center.y
        cz = cz if cz else self.center.z
        sx = sx if sx else self.box_size.x
        sy = sy if sy else self.box_size.y
        sz = sz if sz else self.box_size.z
        content = f'''
receptor = {receptor_path}
ligand = {ligand_path}

center_x = {cx}
center_y = {cy}
center_z = {cz}

size_x = {sx}
size_y = {sy}
size_z = {sz}


exhaustiveness = {exhaustiveness}

num_modes = {num_modes}

energy_range = {energy_range}

out = {outpdbqtfile_path}
'''
        file, file_dir = Path(file), Path(file).parent
        if not file_dir.exists(): file_dir.mkdir(parents=True)
        with open(file.absolute().as_posix(), 'w', encoding='utf-8') as f:
            f.write(content)