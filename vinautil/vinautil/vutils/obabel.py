from dataclasses import dataclass, field
from pathlib import Path
from openbabel import pybel
from typing import Optional, NoReturn, List, Union, Dict, Tuple
import re
import tempfile, os
import numpy as np


@dataclass
class PDBQTparser:
    file: Path
    module_lst: List[List[str]] = field(init=False)
    fmt: str = field(default_factory=lambda: 'pdbqt')

    def __post_init__(self):
        self.module_lst = self.get_modules()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.file})"

    def __str__(self):
        return f"offer some method for {self.file}"

    def get_modules(self, pose_num: int = None) -> Union[List[str], str]:
        module_lst = []
        module = []
        txt = self.file.read_text()
        module_lst = re.findall(r'MODEL\s+\d+.*?ENDMDL', txt, re.DOTALL)
        if pose_num is None:
            return module_lst
        else:
            return module_lst[pose_num]

    @property
    def module_num(self) -> int:
        return len(self.module_lst)

    @property
    def scores(self) -> Dict[str,float]:# get affinity from pdbqt file
        vina_match = r"VINA RESULT:\s+(-\d+\.\d+)"
        model_match = r"MODEL\s+(\d+)"
        return {str(re.search(model_match, i).group(1)): float(re.search(vina_match, i).group(1)) for i in self.module_lst}

    @property
    def score_density(self) -> Dict[str,float]:
        return {k: v/self.atomnum() for k,v in self.scores.items()}

    @property
    def rerank_scores(self) -> List[Tuple[float, int]]:
        l = self.scores.values()
        s = sorted(list(set(l)), reverse=False)
        return [(i, ii) for i in l for ii,jj in enumerate(s) if i == jj]

    @property
    def energies(self):
        return [float(i) for i in self.scores.values()]

    @staticmethod
    def get_coord_dict(fmt: str, file: Path, removeHs: bool = True) -> dict:
        molH = next(pybel.readfile(format=fmt, filename=file.as_posix())) # read first molecule
        if removeHs: molH.OBMol.DeleteHydrogens()
        return {atom.idx: atom.coords for atom in molH}

    def atomnum(self, removeHs: bool = True) -> int:
        mols = pybel.readfile(format=self.fmt, filename=self.file.as_posix())
        if removeHs:
            mols = [i.OBMol.DeleteHydrogens() and i for i in mols]
        else:
            mols = [i for i in mols]
        atom_counts = [len(mol.atoms) for mol in mols]
        if all(count == atom_counts[0] for count in atom_counts):
            # print("All structures have the same number of atoms.")
            return atom_counts[0]
        else:
            raise ValueError("Structures have different numbers of atoms.")
      
    @staticmethod
    def pose_coords(pdbqt_string: str, remove_H = True) -> np.array:
        mol_pose = pybel.readstring('pdbqt', pdbqt_string)
        if remove_H: mol_pose.OBMol.DeleteHydrogens()
        # 创建一个空的numpy数组，用来存储坐标信息
        coords = np.empty((len(mol_pose.atoms), 3))
        # 遍历mol_pose1的原子，把每个原子的坐标赋值给coords数组
        for i, atom in enumerate(mol_pose.atoms):
            coords[i] = atom.coords
        return coords
    
    @staticmethod
    def rmsd(pose_ref:np.array, pose:np.array):
        # 计算两个数组之间的差值，然后对每个元素平方
        diff = pose_ref - pose
        diff_squared = diff ** 2
        assert len(pose_ref) == len(pose)
        N = len(pose_ref)
        # 对平方后的数组求和，然后除以原子个数N，最后开平方根
        sum_diff_squared = diff_squared.sum ()
        rmsd = np.sqrt (sum_diff_squared / N)
        return rmsd



@dataclass
class PDBQTtoMol2:
    '''
    pdbqt 还原为mol2格式，还原健价信息, 删除所有氢元素方便后续计算RMSDß
    '''
    original_mol2_file: str
    undock_pdbqt: str
    docked_pdbqt: List[str] # autodock vina outputs
    pybel_mol: (List[pybel.Molecule], type(None)) = field(init=False)

    def __post_init__(self):
        docked_pdbqt_pybel_mol = [pybel.readstring('pdbqt', i) for i in self.docked_pdbqt]
        self.pybel_mol = self.get_pybel_mol(docked_pdbqt_pybel_mol)
    
    @staticmethod
    def get_coord_dict(fmt: str, file_string: str, removeHs: bool = True) -> dict:
        molH = pybel.readstring(fmt, file_string) # read first molecule
        if removeHs: molH.OBMol.DeleteHydrogens()
        return {atom.idx: atom.coords for atom in molH}
    
    def get_pybel_mol(self, docked_pdbqt_pybel_mol) -> List[pybel.Molecule]:
        '''
        mind refence https://github.com/ag83/pdbqt-to-mol2/blob/be40bdda20ffb96cd3d173accf77e7a2da9a49aa/convert_to_mol2.py#L15
        convert autodock vina dock results restore to mol2 format, which include bond infomations
        :return: pybel.Molecule or None
        '''
        undocked_pdbqt = self.get_coord_dict('pdbqt', self.undock_pdbqt)
        original_mol2 = self.get_coord_dict('mol2', self.original_mol2_file)
        if len(docked_pdbqt_pybel_mol) == 1:
            i_coords = {atom.idx: atom.coords for atom in docked_pdbqt_pybel_mol[0]}
            mol = self.__update_coordinates(docked_pdbqt=i_coords, undocked_pdbqt=undocked_pdbqt, original_mol2=original_mol2)
            mol.OBMol.DeleteHydrogens()
            return [mol]
        elif len(docked_pdbqt_pybel_mol) > 1:
            lst = []
            for i in docked_pdbqt_pybel_mol:
                i.OBMol.DeleteHydrogens()
                i_coords = {atom.idx: atom.coords for atom in i}
                lst.append(self.__update_coordinates(docked_pdbqt=i_coords, undocked_pdbqt=undocked_pdbqt, original_mol2=original_mol2))
            return lst
        else:
            return []

    def __update_coordinates(self, undocked_pdbqt: Dict, docked_pdbqt: Dict, original_mol2: Dict) -> Optional[pybel.Molecule]:
        assert len(undocked_pdbqt) == len(docked_pdbqt) == len(
            original_mol2), f'Not equal number of atoms in molecules\n' \
                            f'undocked_pdbqt: {len(undocked_pdbqt)}\n' \
                            f'docked_pdbqt: {len(docked_pdbqt)}\n' \
                            f'original_mol2: {len(original_mol2)}'
        original_coord = {}
        for key in original_mol2:
            coord_update = [round(x, 3) for x in original_mol2[key]]
            coord_update = tuple(coord_update)
            original_coord.update({key: coord_update})
        coord_map = {}
        for idx, coord in original_coord.items():
            # potential bottleneck for large molecules
            for ind, coordinates in undocked_pdbqt.items():
                n = 0
                if coord[0] == coordinates[0]:
                    n = n + 1
                    if coord[1] == coordinates[1]:
                        n = n + 1
                        if coord[2] == coordinates[2]:
                            n = n + 1
                    else:
                        if coord[2] == coordinates[2]:
                            n = n + 1
                else:
                    if coord[1] == coordinates[1]:
                        n = n + 1
                        if coord[2] == coordinates[2]:
                            n = n + 1
                    else:
                        if coord[2] == coordinates[2]:
                            n = n + 1
                if n == 3:
                    coord_map.update({idx: ind})
                elif n == 2:
                    if idx in coord_map:
                        pass
                    else:
                        coord_map.update({idx: ind})
                elif n == 1:
                    if idx in coord_map:
                        pass
                    else:
                        coord_map.update({idx: ind})
                else:
                    pass
        if len(coord_map) == len(original_mol2):
            coord_conform = {}
            for index1, index2 in coord_map.items():
                coord_conform.update({index1: docked_pdbqt.get(index2)})
            mol2 = pybel.readstring('mol2', self.original_mol2_file)
            mol2.OBMol.DeleteHydrogens()
            for atom in mol2:
                atom.OBAtom.SetVector(coord_conform.get(atom.idx)[0], coord_conform.get(atom.idx)[1],
                                      coord_conform.get(atom.idx)[2])
            self.pybel_mol = mol2
            return mol2
        else:
            raise ValueError('Lost coordinates in mapping')

    def to_file(self, out_file: Path, fmt: str = 'mol2') -> NoReturn:
        if isinstance(self.pybel_mol, list):
            for i in self.pybel_mol:
                with pybel.Outputfile(fmt, out_file.as_posix(), overwrite=True) as f:
                    f.write(i)
        elif isinstance(self.pybel_mol, pybel.Molecule):
            self.pybel_mol.write(fmt, out_file.as_posix(), overwrite=True)
        else:
            raise print('No pybel.Molecule found')

    def to_string(self, fmt: str = 'mol2') -> Union[List[str], str]:
        if isinstance(self.pybel_mol, list):
            return [i.write(fmt) for i in self.pybel_mol]
        elif isinstance(self.pybel_mol, pybel.Molecule):
            return self.pybel_mol.write(fmt)
        else:
            raise print('No pybel.Molecule found')

def mol2pdb(read_file: Path, out_file: Path, s_fmt='mol2', t_fmt='pdb' ):
    # mol2 文件格式转化为pdb格式
    mols = list(pybel.readfile(format=s_fmt, filename=read_file.as_posix()))
    if len(mols) == 1:
        mol = mols[0]
        mol.write(t_fmt, out_file.as_posix(), overwrite=True)

