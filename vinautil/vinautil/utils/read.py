#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :read
@Description:       : read cheminfo file
@Date               :2023/4/9 11:14:29
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''

import re
import pandas as pd
from utils.pdbparse_module import Bdp

def read_mol2_file(file_path):
    # 读取文件内容
    with open(file_path, 'r') as file:
        content = file.read()
    # 使用正则表达式匹配@<TRIPOS>ATOM和下一个@之间的内容
    atom_section = re.search(r'@<TRIPOS>ATOM(.*?)@', content, re.DOTALL)
    # 如果找到匹配项，提取原子记录并将其分割成单独的行
    if atom_section:
        atom_records = atom_section.group(1).strip().split('\n')
    else:
        raise ValueError('未找到@<TRIPOS>ATOM节')
    # 将原子记录分割成单独的字段，并将其保存为pandas DataFrame
    atom_data = [record.split() for record in atom_records]
    df = pd.DataFrame(atom_data, columns=['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge'])
    return df

def get_mol_covalent_atom_coords(df:pd.DataFrame, search_atom_name:str): # locate the atom coords in the mol2 file
    df_c = df[df['atom_name'] == search_atom_name]
    if len(df_c.index) == 1:
        ind = df_c.index[0]
        x,y,z = df_c.loc[ind, 'x'], df_c.loc[ind, 'y'], df_c.loc[ind, 'z']
        return x,y,z

def get_rec_covalent_atom_coords(pdbfile, resname: str, chain: str, resseq: str, atomname: str):
    # locate the atom coords in the pdb file
    df = Bdp.clean_struct(file=pdbfile, label='ATOM')
    # 去除所有字符串列中的前后空格
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    bdp = df
    # keep locate the atom coords failed in the pdb file, clean alpha construct in the pdb file
    df_c = bdp[(bdp['ResidueName'] == resname.upper()) & (bdp['ChainIndentifier'] == chain.upper()) & (bdp['ResidueSequence'] == str(resseq)) & (bdp['AtomName'] == atomname)]
    if len(df_c.index) == 1:
        ind = df_c.index[0]
        x,y,z = df_c.loc[ind, 'X'], df_c.loc[ind, 'Y'], df_c.loc[ind, 'Z']
        return x,y,z


