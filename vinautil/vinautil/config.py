#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :configutils.py
@Description:       : generate config file for autodock vina
@Date     :2022/12/28 21:58:08
@Author      :hotwa
@version      :1.1
'''
import sys
from pathlib import Path
here = Path(__file__).parent.resolve()
from dataclasses import dataclass, field
from typing import Iterable
from openbabel import pybel
try:
    from .typecheck import typeassert
except ImportError:
    from utils.typecheck import typeassert

@typeassert(file=Path, fmt=str)
class molecule:
    file: Path

    def __init__(self,file: Path, fmt:str=None):
        self.file = file
        self.fmt = fmt if bool(fmt) else file.suffix.replace('.','')

    def get_coord_lst(self, addHs:bool = True):
        molH = pybel.readfile(format = self.fmt, filename=self.file.__str__())
        molH = next(molH)
        # molH.OBMol.DeleteHydrogens()
        if addHs: molH.OBMol.AddPolarHydrogens()
        return [list(atom.coords) for atom in molH]

    def get_coord_dict(self, addHs:bool = True):
        molH = pybel.readfile(format = self.fmt, filename=self.file.__str__())
        molH = next(molH)
        if addHs: molH.OBMol.AddPolarHydrogens()
        return {atom.idx: atom.coords for atom in molH}

    def compute_box(self, extending = 10):
        coords_lst = self.get_coord_lst()
        xlst, ylst, zlst = [i[0] for i in coords_lst], [i[1] for i in coords_lst], [i[2] for i in coords_lst]
        ([minX, minY, minZ], [maxX, maxY, maxZ]) = ([min(xlst), min(ylst), min(zlst)], [max(xlst), max(ylst), max(zlst)])
        minX = minX - float(extending)
        minY = minY - float(extending)
        minZ = minZ - float(extending)
        maxX = maxX + float(extending)
        maxY = maxY + float(extending)
        maxZ = maxZ + float(extending)
        SizeX = maxX - minX
        SizeY = maxY - minY
        SizeZ = maxZ - minZ
        CenterX = (maxX + minX) / 2
        CenterY = (maxY + minY) / 2
        CenterZ = (maxZ + minZ) / 2
        c = {
                'x': float("{:.2f}".format(CenterX)),
                'y': float("{:.2f}".format(CenterY)),
                'z': float("{:.2f}".format(CenterZ)),
        }
        s = {
                'x': float("{:.2f}".format(SizeX)),
                'y': float("{:.2f}".format(SizeY)),
                'z': float("{:.2f}".format(SizeZ))
        }
        rd = {
            'center': xyz_point(**c),
            'box_size': xyz_point(**s)
        }
        return rd

@typeassert(x=(float, int), y=(float, int), z=(float, int))
@dataclass
class xyz_point:
    x: (float, int)
    y: (float, int)
    z: (float, int)

    def __repr__(self):
        return f"xyz_point(x={self.x}, y={self.y}, z={self.z})"

    def __str__(self):
        return f"three demensional point: x={self.x}, y={self.y}, z={self.z}"

    @classmethod
    def from_iterable(cls,iterable_object:Iterable):
        if isinstance(iterable_object,Iterable):
            return cls(*iterable_object)
        else:
            raise TypeError('iterable_object must be Iterable')

@typeassert(receptor=Path,
    ligand=Path,
    box_size=(xyz_point, type(None)),
    center=(xyz_point, type(None)),
    outpdbqtfile=Path,
    exhaustiveness=int,
    num_modes=int,
    energy_range=int)
@dataclass
class vinaconfig:
    receptor: Path
    ligand: Path
    outpdbqtfile: Path
    center: (xyz_point, type(None)) = field(default_factory=lambda: None)
    box_size: (xyz_point, type(None)) = field(default_factory=lambda: None)
    exhaustiveness: int = field(default_factory=lambda: 32)
    num_modes: int = field(default_factory=lambda: 20)
    energy_range: int = field(default_factory=lambda: 5)
    extend: int = field(default_factory=lambda: 10)

    def __str__(self):
        return f"""the class of autodock vina config file:\n\nreceptor: {self.receptor}\nligand: {self.ligand}\ncenter: {self.center}\nbox_size: {self.box_size}\noutpdbqtfile: {self.outpdbqtfile}\nexhaustiveness: {self.exhaustiveness}\nnum_modes: {self.num_modes}\nenergy_range: {self.energy_range}\nextend: {self.extend}\n"""

    def __post_init__(self):
        if (self.center == None) and (self.box_size == None):
            m = molecule(self.receptor)
            box = m.compute_box(extending = self.extend)
            self.center = box['center']
            self.box_size = box['box_size']

    def to_dict(self)->dict:
        return {'receptor': self.receptor,
                'ligand': self.ligand,
                'center': self.center,
                'box_size': self.box_size,
                'outpdbqtfile': self.outpdbqtfile,
                'exhaustiveness': self.exhaustiveness,
                'num_modes': self.num_modes,
                'energy_range': self.energy_range}

    def to_txt(self, file:Path,
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
        file,file_dir = Path(file), Path(file).parent
        if not file_dir.exists(): file_dir.mkdir(parents=True)
        with open(file.absolute().as_posix(), 'w', encoding='utf-8') as f:
            f.write(content)