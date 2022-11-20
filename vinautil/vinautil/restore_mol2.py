#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        : restore_mol2.py
@Description:       : restore docked pdbqt to mol2
@Date     :2022/11/17 10:26:36
@Author      :hotwa
@version      :1.0
'''
import os
from pathlib import Path, PurePath
from openbabel import pybel

def get_coord_dict(fmt, file):
    molH = [i for i in pybel.readfile(format = fmt, filename=file.__str__())][0]
    molH.OBMol.DeleteHydrogens()
    return {atom.idx: atom.coords for atom in molH}

def pdbqt2mol2(original_mol2_file, undock_pdbqt, docked_pdbqt, out_mol2):
    # mind refence https://github.com/ag83/pdbqt-to-mol2/blob/be40bdda20ffb96cd3d173accf77e7a2da9a49aa/convert_to_mol2.py#L15
    # convert autodock vina dock results restore to mol2 format, which include bond infomations
    undocked_pdbqt = get_coord_dict('pdbqt', undock_pdbqt)
    docked_pdbqt = get_coord_dict('pdbqt', docked_pdbqt)
    original_mol2 = get_coord_dict('mol2', original_mol2_file)
    assert len(undocked_pdbqt) == len(docked_pdbqt) == len(original_mol2),f'Not equal number of atoms in molecules\n{docked_pdbqt}'
    original_coord = {}
    for key in original_mol2:
        coord_update = [round(x, 3) for x in original_mol2[key]]
        coord_update = tuple(coord_update)
        original_coord.update({key: coord_update})
    coord_map = {}
    for idx, coord in original_coord.items():
        # potential bottleneck for large molecules
        for ind, coordinates in undocked_pdbqt.items():
            n=0
            if coord[0] == coordinates[0]:
                n=n+1
                if coord[1] == coordinates[1]:
                    n=n+1
                    if coord[2] == coordinates[2]:
                        n=n+1
                else:
                    if coord[2] == coordinates[2]:
                        n=n+1
            else:
                if coord[1] == coordinates[1]:
                    n=n+1
                    if coord[2] == coordinates[2]:
                        n=n+1
                else:
                    if coord[2] == coordinates[2]:
                        n=n+1
            if n == 3:
                coord_map.update({idx:ind})
            elif n == 2:
                if idx in coord_map:
                    pass
                else:
                    coord_map.update({idx:ind})
            elif n == 1:
                if idx in coord_map:
                    pass
                else:
                    coord_map.update({idx:ind})
            else:
                pass
    if len(coord_map) == len(original_mol2):
        coord_conform = {}
        for index1, index2 in coord_map.items():
            coord_conform.update({index1:docked_pdbqt.get(index2)})
        mol2 = pybel.readfile('mol2', original_mol2_file)
        mol2 = next(mol2)
        mol2.OBMol.DeleteHydrogens()
        for atom in mol2:
            atom.OBAtom.SetVector(coord_conform.get(atom.idx)[0], coord_conform.get(atom.idx)[1], coord_conform.get(atom.idx)[2])
        mol2.write('mol2', out_mol2, overwrite=True)
    else:
        print('Lost coordinates in mapping')

if __name__ == '__main__':
    # restore docked_pdbqt to mol2
    docked_pdbqt = f'test_file/test_mgltools_vina112_dock_2xd1/vina112_results/2xd1_A-2XD1_A_CEF/2xd1_A-2XD1_A_CEF_0.pdbqt'
    undocked_pdbqt = f'test_file/test_mgltools_vina112_dock_2xd1/pymol_addHs_mgltools_CEF.pdbqt'
    original_mol2 = f'test_file/pymol_addHs_CEF.mol2'
    out_mol2 = f'test_file/restore_CEF.mol2'
    pdbqt2mol2(original_mol2, undocked_pdbqt, docked_pdbqt, out_mol2)