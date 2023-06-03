#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :ledock_utils.py
@Description:       : ledock准备工作,分割任务
@Date               :2023/2/22 18:47:25
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
from pathlib import Path
from typing import List
import subprocess

def generate_ledock_file(receptor:Path, l_list_outfile:Path, out_file:Path, l_list: Path,
                         x=[0, 0], y=[0, 0], z=[0, 0], n_poses=20,rmsd=1.0):
    rmsd = str(rmsd)
    x = [str(x) for x in x]
    y = [str(y) for y in y]
    z = [str(z) for z in z]
    n_poses = str(n_poses)
    l_list_outfile.write_text(l_list.as_posix()) # only for one molecular docking, otherwise change this part

    file = [
        'Receptor\n',
        receptor.as_posix() + '\n\n',
        'RMSD\n',
        rmsd + '\n\n',
        'Binding pocket\n',
        x[0], ' ', x[1], '\n',
        y[0], ' ', y[1], '\n',
        z[0], ' ', z[1], '\n\n',
        'Number of binding poses\n',
        n_poses + '\n\n',
        'Ligands list\n',
        l_list_outfile.as_posix() + '\n\n',
        'END']
    print(f'save ledock config: {out_file.as_posix()}')
    out_file.write_text(''.join(file))

def box_vina_to_ledock(center_x, center_y, center_z, size_x, size_y, size_z):
    # 计算 LeDockBox 中的坐标值
    minX = round(center_x - size_x / 2, 2)
    maxX = round(center_x + size_x / 2, 2)
    minY = round(center_y - size_y / 2, 2)
    maxY = round(center_y + size_y / 2, 2)
    minZ = round(center_z - size_z / 2, 2)
    maxZ = round(center_z + size_z / 2, 2)

    # 构造 LeDockBox 字符串
    LeDockBox = "*********LeDock Binding Pocket*********\n" + \
                "Binding pocket\n%.1f %.1f\n%.1f %.1f\n%.1f %.1f\n" % (minX, maxX, minY, maxY, minZ, maxZ)
    print(LeDockBox)

    return minX, maxX, minY, maxY, minZ, maxZ

def config_vina2ledock(vina_file: Path, out_path: Path):
    # test one instance
    if not vina_file.exists(): raise FileNotFoundError("File not found: " + vina_file.as_posix())
    # extract box infos
    content_list = [i for i in vina_file.read_text().splitlines() if i]  # read vina config
    center_x, center_y, center_z, size_x, size_y, size_z = [float(content.split(" = ")[1]) for content in
                                                            content_list[2:8]]
    center_x = round(center_x, 2)
    center_y = round(center_y, 2)
    center_z = round(center_z, 2)
    size_x = round(size_x, 2)
    size_y = round(size_y, 2)
    size_z = round(size_z, 2)
    minX, maxX, minY, maxY, minZ, maxZ = box_vina_to_ledock(center_x, center_y, center_z, size_x, size_y, size_z)
    # extract ligand infos
    # 去掉_to_pdbqt 将_ploar.pdbqt变为_polar.pdb
    lig_content = content_list[1].split('=')[1].strip().replace('molecular_polar_mol2_to_pdbqt',
                                                                'molecular_polar_mol2', ).replace('_ploar.pdbqt',
                                                                                                  '_polar.mol2')
    lig_content = out_path.parent.joinpath(lig_content)
    if not lig_content.exists(): raise FileNotFoundError("File not found: " + lig_content.as_posix())
    # extract protein infos
    pro_content = out_path.parent.joinpath(content_list[0].split('=')[1].strip().replace('_pdbqt', '').replace('pdbqt', 'pdb'))
    if not pro_content.exists(): raise FileNotFoundError("File not found: " + pro_content.as_posix())
    liglst = out_path.joinpath(vina_file.stem + '.list')
    ledock_config_file = out_path.joinpath(vina_file.stem + '.in')
    generate_ledock_file(receptor=pro_content, x=[minX, maxX], y=[minY, maxY], z=[minZ, maxZ], n_poses=20,
                         l_list=lig_content, l_list_outfile=liglst,
                         out_file=ledock_config_file)

def ledock_split(ledock_results: Path, option='-spli'):
    ledock_binary = '/home/lyzeng/ledock_benchmark/ledock_linux_x86'
    res = subprocess.run([ledock_binary, option, ledock_results.as_posix()])
    if res.returncode != 1:
        return res

def split_main():
    ledock_results = list(Path('/home/lyzeng/ledock_benchmark/molecular_polar_mol2').glob('*.dok'))
    file = ledock_results[0]
    split_res = list(map(lambda x: x!=None ,[ledock_split(i) for i in ledock_results]))

if __name__ == '__main__':
    ...




