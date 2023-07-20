#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :deeppocket
@Description:       :
@Date               :2023/7/15 11:09:36
@Author             :lyzeng
@mail               :pylyzeng@gmail.com
@version            :1.0
'''
from dataclasses import dataclass, field
import pandas as pd
import deepchem as dc
import mdtraj as md
from pymol import cmd
from pathlib import Path
import nglview as nv
from time import sleep
import py3Dmol
import numpy as np
from typing import List, Tuple
from pdbfixer import PDBFixer
from openmm.app import PDBFile

@dataclass
class AtomInfo:
    name: str
    resn: str
    resi: int
    chain: str
    coord: List[float]

@dataclass
class deeppocket:
    file: Path
    pid: str = field(init=False)

    def __post_init__(self):
        self.pid = self.file.stem

    @property
    def PocketDataFrame(self):
        if not hasattr(self, '_PocketDataFrame'):
            self._PocketDataFrame = self.get_pocket()
        return self._PocketDataFrame

    @property
    def file_clean(self):
        if not hasattr(self, '_file_clean'):
            self._file_clean = self.cleanATOM()
        return self._file_clean

    @property
    def file_pdbfixer(self):
        if not hasattr(self, '_file_pdbfixer'):
            self._file_pdbfixer = self.fix_pdb()
        return self._file_pdbfixer

    @property
    def atoms(self):
        if not hasattr(self, '_atoms'):
            self._atoms = self.get_atoms(file=self.file)
        return self._atoms

    @property
    def chains(self):
        if not hasattr(self, '_chains'):
            self._chains = self.get_chains()
        return self._chains

    @property
    def hetatm(self):
        if not hasattr(self, '_hetatm'):
            self._hetatm = self.get_hetatm()
        return self._hetatm

    def __str__(self):
        return f"use deepchem search pocket for {self.pid}"

    def get_chains(self) -> list:
        cmd.reinitialize('everything')
        cmd.load(self.file.as_posix())
        return cmd.get_chains(selection=self.pid)

    def get_hetatm(self) -> list:
        cmd.reinitialize('everything')
        cmd.load(self.file.as_posix())
        return self.pymol_get_resn('hetatm')

    def get_atoms(self, file: Path, selection='all'):
        cmd.reinitialize('everything')
        cmd.load(file.as_posix())
        model = cmd.get_model(selection)
        atoms = []
        for atom in model.atom:
            atom_info = AtomInfo(
                name=atom.name,
                resn=atom.resn,
                resi=atom.resi,
                chain=atom.chain,
                coord=[atom.coord[0], atom.coord[1], atom.coord[2]]
            )
            atoms.append(atom_info)
        cmd.reinitialize('everything')
        return atoms

    def cleanATOM(self, out_file=None, ext="_clean.pdb") -> Path:  # from pyrosetta.toolbox import cleanATOM
        """Extract all ATOM and TER records in a PDB file and write them to a new file.

        Args:
            pdb_file (str): Path of the PDB file from which ATOM and TER records
                will be extracted
            out_file (str): Optional argument to specify a particular output filename.
                Defaults to <pdb_file>.clean.pdb.
            ext (str): File extension to use for output file. Defaults to ".clean.pdb"
        """
        # find all ATOM and TER lines
        with open(self.file.as_posix(), "r") as fid:
            good = [l for l in fid if l.startswith(("ATOM", "TER"))]

        # default output file to <pdb_file>_clean.pdb
        if out_file is None:
            out_file = self.file.stem + ext

        # write the selected records to a new file
        with open(out_file, "w") as fid:
            fid.writelines(good)
        return Path(out_file)

    def get_atoms_in_box(self, row: pd.core.series.Series) -> List[AtomInfo]:
        # 从行中获取盒子的范围
        x_start, x_end = row['x_start'], row['x_end']
        y_start, y_end = row['y_start'], row['y_end']
        z_start, z_end = row['z_start'], row['z_end']

        # 筛选出在盒子范围内的原子
        atoms_in_box = [atom for atom in self.atoms
                        if x_start <= atom.coord[0] <= x_end
                        and y_start <= atom.coord[1] <= y_end
                        and z_start <= atom.coord[2] <= z_end]

        return atoms_in_box

    def get_pocket_atoms(self, pocket_residues):
        return list(filter(lambda atom: atom.resi in pocket_residues, self.atoms))

    def get_atom_info(self, atom_name):
        return [atom for atom in self.atoms if atom.name == atom_name]

    def pymol_get_resn(self, object_name: str = 'all') -> Tuple[str]:
        resn_list = []
        cmd.iterate(selection=f'{self.pid} & {object_name}', expression=lambda atom: resn_list.append(atom.resn))
        return tuple(set(resn_list))

    def fix_pdb(self, ph: float = 7.0):
        out_file = self.file_clean.parent.joinpath(f'{self.pid}_fixed_pH_{int(ph)}.pdb')
        if not out_file.exists():
            print("Creating PDBFixer...")
            fixer = PDBFixer(self.file_clean.as_posix())
            print("Finding missing residues...")
            fixer.findMissingResidues()

            chains = list(fixer.topology.chains())
            keys = fixer.missingResidues.keys()
            for key in list(keys):
                chain = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    print("ok")
                    del fixer.missingResidues[key]

            print("Finding nonstandard residues...")
            fixer.findNonstandardResidues()
            print("Replacing nonstandard residues...")
            fixer.replaceNonstandardResidues()
            print("Removing heterogens...")
            fixer.removeHeterogens(keepWater=True)
            print("Finding missing atoms...")
            fixer.findMissingAtoms()
            print("Adding missing atoms...")
            fixer.addMissingAtoms()
            print("Adding missing hydrogens...")
            fixer.addMissingHydrogens(7)
            print("Writing PDB file...")
            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                open(out_file, "w"), keepIds=True)
        return Path(out_file)

    def get_pocket(self) -> pd.DataFrame:
        finder = dc.dock.binding_pocket.ConvexHullPocketFinder()
        pockets = finder.find_pockets(self.file.as_posix())

        x_length = lambda x: x.x_range[1] - x.x_range[0]
        y_length = lambda y: y.y_range[1] - y.y_range[0]
        z_length = lambda z: z.z_range[1] - z.z_range[0]

        # 计算每个口袋的中心点坐标
        x_center = lambda x: (x.x_range[1] + x.x_range[0]) / 2
        y_center = lambda y: (y.y_range[1] + y.y_range[0]) / 2
        z_center = lambda z: (z.z_range[1] + z.z_range[0]) / 2

        # 计算每个口袋的尺寸
        x_lengths = [x_length(pocket) for pocket in pockets]
        y_lengths = [y_length(pocket) for pocket in pockets]
        z_lengths = [z_length(pocket) for pocket in pockets]

        # 计算每个口袋的中心点坐标
        x_centers = [x_center(pocket) for pocket in pockets]
        y_centers = [y_center(pocket) for pocket in pockets]
        z_centers = [z_center(pocket) for pocket in pockets]

        x_starts = [pocket.x_range[0] for pocket in pockets]
        x_ends = [pocket.x_range[1] for pocket in pockets]
        y_starts = [pocket.y_range[0] for pocket in pockets]
        y_ends = [pocket.y_range[1] for pocket in pockets]
        z_starts = [pocket.z_range[0] for pocket in pockets]
        z_ends = [pocket.z_range[1] for pocket in pockets]

        # 创建一个包含所有数据的DataFrame
        df = pd.DataFrame({
            'x_length': x_lengths,
            'y_length': y_lengths,
            'z_length': z_lengths,
            'x_center': x_centers,
            'y_center': y_centers,
            'z_center': z_centers,
            'x_start': x_starts,
            'x_end': x_ends,
            'y_start': y_starts,
            'y_end': y_ends,
            'z_start': z_starts,
            'z_end': z_ends,
            # 'pocket': pockets  # 新增的列，保存Box对象
        })
        df['volume'] = df['x_length'] * df['y_length'] * df['z_length']
        # 按照体积从小到大排序
        df = df.sort_values(by='volume').reset_index(drop=True)
        return df

    def _add_mesh(self, view: nv.widget.NGLWidget, pocket_dict:dict, opacity: float = 0.2, color_scheme: str = 'electrostatic'):
        # view.add_surface(selection='residue number range: 3-40', color='residueindex', opacity=0.2)
        # view.add_surface(selection='chainname:A and 10', color_scheme="electrostatic", opacity=0.2)
        if len(pocket_dict.keys()) != 1: # 口袋有多个链
            print(f'Warning: more than one chain in pocket, chain id: {pocket_dict.keys()}')
        for chain in pocket_dict.keys():
            if chain in self.chains:
                pocket_residue_number = 'or'.join(pocket_dict[chain])
                selection = f":{chain} and ({pocket_residue_number})"
                view.add_surface(selection=selection, color_scheme=color_scheme, opacity=opacity)
        return view



    def draw_box_nglview(self, row: int, color: list = [0, 1, 1], cylinder_radius: float = 0.1, add_mesh: bool = True,
                         opacity: float = 0.3, color_scheme: str = 'electrostatic'):
        """
        draw a box for a pocket
        Args:
            row: line number of the pocket dataframe
            color: a list of three numbers representing the color of the box in RGB format
            opacity: a number representing the opacity of the box

        Returns:
            nglview object
        """
        # 使用MDTraj加载PDB文件
        traj = md.load(self.file.as_posix())

        # 使用nglview创建一个视图对象
        view = nv.show_mdtraj(traj)

        # 从df中获取口袋的信息
        pocket = self.PocketDataFrame.iloc[row]

        # 计算口袋的各个角的坐标
        x_start, y_start, z_start = pocket['x_start'], pocket['y_start'], pocket['z_start']
        x_end, y_end, z_end = pocket['x_end'], pocket['y_end'], pocket['z_end']

        # 定义长方体的16个边的端点
        edges = [
            [[x_start, y_start, z_start], [x_start, y_start, z_end]],
            [[x_start, y_start, z_start], [x_start, y_end, z_start]],
            [[x_start, y_start, z_start], [x_end, y_start, z_start]],
            [[x_start, y_end, z_end], [x_start, y_start, z_end]],
            [[x_start, y_end, z_end], [x_start, y_end, z_start]],
            [[x_start, y_end, z_end], [x_end, y_end, z_end]],
            [[x_end, y_start, z_end], [x_start, y_start, z_end]],
            [[x_end, y_start, z_end], [x_end, y_start, z_start]],
            [[x_end, y_start, z_end], [x_end, y_end, z_end]],
            [[x_end, y_end, z_start], [x_start, y_end, z_start]],
            [[x_end, y_end, z_start], [x_end, y_start, z_start]],
            [[x_end, y_end, z_start], [x_end, y_end, z_end]],
            [[x_start, y_end, z_start], [x_start, y_start, z_start]],
            [[x_start, y_end, z_start], [x_end, y_end, z_start]],
            [[x_end, y_end, z_end], [x_start, y_end, z_end]],
            [[x_end, y_end, z_end], [x_end, y_start, z_end]]
        ]

        # 添加圆柱体来表示长方体的边
        for start, end in edges:
            view.shape.add_cylinder(start, end, color, cylinder_radius)

        if add_mesh:
            atoms = self.get_atoms_in_box(pocket)
            chains = list(set([i.chain for i in atoms]))
            pocket_dict = {chain: [atom.resi for atom in atoms if atom.chain == chain] for chain in chains}
            self._add_mesh(view, pocket_dict, opacity=opacity, color_scheme=color_scheme)

        return view

    def draw_box_py3Dmol(self, row: int, color: str = 'red', opacity: float = 0.25) -> py3Dmol.view:
        """
        draw a box for a pocket using py3Dmol
        Args:
            row: line number of the pocket dataframe

        Returns:
            py3Dmol view object
        """
        df = self.PocketDataFrame
        pocket = df.loc[row]
        # 首先，加载PDB文件
        with open(self.file, 'r') as f:
            pdb_data = f.read()
        # 然后，初始化一个Py3Dmol view
        p = py3Dmol.view()
        # 将PDB数据添加到视图中
        p.addModel(pdb_data, 'pdb')
        # 添加cartoon样式显示蛋白质
        p.setStyle({'cartoon': {'color': 'spectrum'}})
        p.addBox({'center': {'x': pocket.x_center, 'y': pocket.y_center, 'z': pocket.z_center},
                  'dimensions': {'w': pocket.x_length, 'h': pocket.y_length, 'd': pocket.z_length},
                  'color': color, 'opacity': opacity})

        # 显示分子和口袋
        p.zoomTo()
        return p
