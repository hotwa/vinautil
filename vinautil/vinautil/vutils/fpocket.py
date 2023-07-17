#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :fpocket
@Description:       :
@Date               :2023/7/17 16:15:29
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
from dataclasses import dataclass, field
import math
import re, os
from pathlib import Path
import numpy as np
import mdtraj as md
import nglview as nv
import subprocess
import datetime
from shutil import rmtree


@dataclass
class fpocket_calulate:
    pdb_file: Path
    overwrite: bool = field(default=True)
    conda_prefix: Path = field(default_factory=lambda: Path(os.environ.get('CONDA_PREFIX')))
    fpocket: Path = field(default_factory=lambda: Path(os.environ.get('CONDA_PREFIX')) / 'bin/fpocket')

    def __post_init__(self):
        self.out_dir = self._search_pocket()
        self.info = self.__get_info()

    @property
    def pocket_number(self):
        return len(self.info)

    def _search_pocket(self):
        fpocket_results = self.pdb_file.parent.joinpath(f'{self.pdb_file.stem}_out')
        if not fpocket_results.exists:
            self._subprocess_fpocket()
            return fpocket_results
        else:
            print(f'{fpocket_results.as_posix()} already exists!')
            if self.overwrite:
                print(f'remove {fpocket_results.as_posix()}')
                rmtree(fpocket_results.as_posix())
                self._subprocess_fpocket()
                return fpocket_results

    def _subprocess_fpocket(self):
        CMD_ = f'{self.fpocket.as_posix()} -f {self.pdb_file.as_posix()}'
        p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
        while p.poll() is None:  # progress still runing
            subprocess_read_res = p.stdout.read().decode('utf-8')
            print(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
        print(f'Analyse {self.pdb_file.as_posix()} sucess!')

    def __get_info(self):
        info_txt = list(self.out_dir.glob('*_info.txt'))[0]
        # 读取_info.txt文件
        with open(info_txt.as_posix(), "r") as file:
            content = file.read()
        # 使用正则表达式匹配口袋信息并提取数值
        pocket_matches = re.findall(r"Pocket (\d+) :\s+(.*?)\n\n", content, re.DOTALL)
        # 定义一个空列表用于存储口袋信息
        pockets = []
        # 遍历口袋匹配结果
        for match in pocket_matches:
            pocket_index = int(match[0])
            pocket_info = match[1].strip().split("\n")
            # 创建口袋字典
            pocket = {"Index": pocket_index}
            # 提取口袋信息并添加到字典中
            for line in pocket_info:
                key, value = line.split(":")
                pocket[key.strip()] = value.strip()
                if key.strip() == 'Cent. of mass - Alpha Sphere max dist':
                    size = math.ceil(float(value.strip()) * 2)
                    center = self._get_mass_center(pocket_index)
                    pocket['AutoDock Vina Box'] = {'center': center, 'size': [size, size, size]}
            # 将口袋字典添加到口袋列表中
            pockets.append(pocket)
        return pockets

    def _get_mass_center(self, num: int):
        pocket_pdbfile = self.out_dir.joinpath(f'pockets/pocket{num}_atm.pdb')
        with open(pocket_pdbfile, 'r') as file:
            lines = file.readlines()
        # 提取原子坐标
        atoms = []
        for line in lines:
            if line.startswith('ATOM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append([x, y, z])
        # 将坐标转换为numpy数组
        atoms = np.array(atoms)
        # 计算质心
        centroid = np.mean(atoms, axis=0)
        return list(centroid)

    def _get_vina_box(self, num: int):
        return self.info[num - 1]['AutoDock Vina Box']

    def _get_box_view(self, num: int):
        # 使用MDTraj加载分子动力学轨迹
        trajectory = md.load(self.pdb_file.as_posix())

        vina_box = self._get_vina_box(num)
        center = vina_box['center']
        size = vina_box['size']

        x_start = center[0] - size[0] / 2
        y_start = center[1] - size[1] / 2
        z_start = center[2] - size[2] / 2

        x_end = center[0] + size[0] / 2
        y_end = center[1] + size[1] / 2
        z_end = center[2] + size[2] / 2

        # 使用nglview创建一个视图
        view = nv.show_mdtraj(trajectory)

        # 显示蛋白质
        view.add_representation('cartoon', selection='protein')

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
        color = [0, 1, 1]
        cylinder_radius = 0.1
        for start, end in edges:
            view.shape.add_cylinder(start, end, color, cylinder_radius)
        return view

