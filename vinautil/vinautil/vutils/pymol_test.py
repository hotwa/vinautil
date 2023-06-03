#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from pathlib import Path
from pymol import cmd, finish_launching
import ast
from typing import Optional, Tuple, NoReturn, Dict, List, Iterable, Callable, Any

def get_resn(name:str = 'all'):
    resn_list = []
    cmd.iterate(selection=name, expression=lambda atom: resn_list.append(atom.resn))
    return tuple(set(resn_list))

def pymol_print(object_name: str = 'hetatm', atom_attribute: str = 'resn') -> Tuple[int, List[str]]:
    """
        打印 pymol 中对象的属性。

        这是一个用来打印 pymol 中对象的属性的函数。

        :param object_name: 对象的名称，默认为 'hetatm'。
        :type object_name: str，optional
        :param atom_attribute: 属性的名称，默认为 'resn'。
        Optional: ['ID', 'alt', 'b', 'cartoon', 'chain', 'color',
        'elec_radius', 'elem', 'flags', 'formal_charge',
        'geom', 'index', 'label', 'model', 'name', 'numeric_type',
        'oneletter', 'p', 'partial_charge', 'protons', 'q', 'rank', 'reps', 'resi',
        'resn', 'resv', 's', 'segi', 'ss', 'stereo', 'text_type', 'type', 'valence', 'vdw']
        :type atom_attribute: str，optional
        :return: 所选对象的数量和属性的值的列表。
    """
    local_namespace = {'res': [], 'cmd': cmd},
    cmd_obj = ast.literal_eval(
        f"cmd.iterate(selection='{object_name}', expression=lambda atom: res.append(atom.{atom_attribute}))"
    )
    select_num = eval(cmd_obj, local_namespace)
    return (select_num, res)

def pymol_print_exec(object_name: str = 'hetatm', atom_attribute: str = 'resn') -> Tuple[int, List[str]]:
    """
        打印 pymol 中对象的属性。

        这是一个用来打印 pymol 中对象的属性的函数。

        :param object_name: 对象的名称，默认为 'hetatm'。
        :type object_name: str，optional
        :param atom_attribute: 属性的名称，默认为 'resn'。
        Optional: ['ID', 'alt', 'b', 'cartoon', 'chain', 'color',
        'elec_radius', 'elem', 'flags', 'formal_charge',
        'geom', 'index', 'label', 'model', 'name', 'numeric_type',
        'oneletter', 'p', 'partial_charge', 'protons', 'q', 'rank', 'reps', 'resi',
        'resn', 'resv', 's', 'segi', 'ss', 'stereo', 'text_type', 'type', 'valence', 'vdw']
        'p' is not supported. 's' is wrapper.SettingWrapper object
        :type atom_attribute: str，optional
        :return: 所选对象的数量和属性的值的列表。
    """
    local_namespace = {'res': [], 'cmd': cmd}
    CMD = f"cmd.iterate(selection='{object_name}', expression=lambda atom: res.append(atom.{atom_attribute}))"
    exec(CMD, local_namespace)
    return local_namespace['res']

def delete(i_cmd: cmd = None):
    if i_cmd: cmd = i_cmd
    if cmd is i_cmd:
        print('same object')
    cmd.delete("all")

if __name__ == '__main__':
    launch_res = finish_launching(['pymol', '-R']) # default port 9123 localhost
    print(launch_res)
    PDBdata = list(Path('C:\\Users\\Admin\\program\\scarversion\\dock20221224\\PDBdatabase20221231').rglob('*.ent.gz'))
    cmd.load(PDBdata[1234])
    cmd.remove('metal solvent')
    hetatm_lst = get_resn('hetatm')
    organic_lst = get_resn('organic')
    print(f'hetatm: {hetatm_lst}')
    print(f'organic: {organic_lst}')

    attribute_lst = ['ID', 'alt', 'b', 'cartoon', 'chain', 'color', 'elec_radius', 'elem', 'flags', 'formal_charge', 'geom', 'index', 'label', 'model', 'name', 'numeric_type', 'oneletter', 'partial_charge', 'protons', 'q', 'rank', 'reps', 'resi', 'resn', 'resv', 's', 'segi', 'ss', 'stereo', 'text_type', 'type', 'valence', 'vdw']
    res = []
    cmd.iterate(selection='hetatm', expression=lambda atom: res.append(atom.ID))
    print(f"{res[0]}----{len(res)}")
    for attribute in attribute_lst:
        print(attribute)
        res = pymol_print_exec('hetatm', attribute)
        print(f"{res[0]}----{len(res)}")
    res = []
    res = pymol_print_exec('hetatm', '')
    print(res)
    cmd.get_coords('hetatm')
    cmd.set_color()








