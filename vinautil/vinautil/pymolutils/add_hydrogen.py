#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :add_hydrogen.py
@Description:       : add hydrogen in pymol
@Date     :2022/10/24 11:17:11
@Author      :hotwa
@version      :1.0
'''
from pymol import cmd
from openbabel import pybel
"""
conda install -c conda-forge pymol-open-source openbabel --yes
"""

def add_polar_hydrogen(file, out_file,read_fmt = 'mol2', pymol_obj='sele',
    save_format='mol2',tool = 'pymol'):
    """add_polar_hydrogen _summary_

    _extended_summary_

    Arguments:
        file {pathlike} -- input file path
        out_file {pathlike} -- output file

    Keyword Arguments:
        read_fmt {string} -- the format of read file (default: {'mol2'})
        pymol_obj {string } -- add polar hydrogen object name, ligand ID (default: {'sele'})
        save_format {string} -- the format of save file (default: {'mol2'})
        tool {string} -- the tool of add polar hydrogens (default: {'pymol'})

    Returns:
        _type_ -- path or None
    """
    if tool == 'pymol':
        cmd.reinitialize("everything")
        cmd.load(filename=file, format=read_fmt)
        cmd.h_add(pymol_obj) # add all hydrogens in this molecular
        cmd.remove(f"{pymol_obj} & hydro & not nbr. (don.|acc.)") # remove no-polar-hydrogen
        cmd.save(filename=out_file,selection=pymol_obj, format=save_format)
        return out_file
    elif tool == 'pybel':
        omol = list(pybel.readfile(format = 'mol2', filename = file))[0]
        omol.OBMol.AddPolarHydrogens()
        omol.write('mol2',out_file,overwrite=True)
        return out_file
    else:
        return None
