#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :pyrosettafix
@Description:       :
@Date               :2023/3/11 16:54:38
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''

# 使用pyrosetta快速修复蛋白缺失的侧链
from pathlib import Path
import pyrosetta
from pyrosetta import rosetta
from multiprocessing import Pool
from pyrosetta.toolbox import cleanATOM
import argparse, os, sys


def fix_optimize(file: Path, out_file: Path):
    '''
    FastRelax使用快速梯度下降算法，可以在较短时间内对蛋白质进行优化，并且对于结构中的非构象缺陷，例如氢键、离子对、芳香性相互作用和溶剂-蛋白质相互作用等进行优化。
    ref2015是Rosetta程序包中的一个分数函数，它是Rosetta2015中引入的一个新的蛋白质力场，用于蛋白质结构预测和设计。这个力场是从先前的Rosetta力场中提炼出来的，经过了一系列的校正和优化，可以更好地预测蛋白质的折叠构象。在FastRelax中，ref2015可以作为一个可选参数来指定使用哪个力场来进行优化。
    :param file:
    :param out_file:
    :return:
    '''
    # 使用pyrosetta修复蛋白结构
    cleanATOM(file.as_posix())
    file = file.parent.joinpath(f'{file.stem}.clean{file.suffix}')
    # 初始化PyRosetta
    pyrosetta.init()
    # 读入蛋白质结构
    pose = pyrosetta.pose_from_pdb(file.as_posix())
    # fix residue side chain
    scorefxn = pyrosetta.create_score_function('ref2015')
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn)
    relax.apply(pose)
    # 输出修复后的结构
    pose.dump_pdb(out_file.as_posix())

def main_multi(p):
    dire = Path(p)
    if not dire.exists(): raise NotADirectoryError(f'Not found: {p.as_posix()}')
    files = list(dire.rglob('*.pdb'))

    # 创建进程池
    with Pool() as pool:
        # 将fix_optimize()函数应用到每个文件
        results = [pool.apply_async(fix_optimize, args=(file, file.parent.joinpath(file.name.split('.')[0] + '_relaxed.pdb'))) for file in files]
        # 等待所有进程完成
        for result in results:
            result.wait()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optimize pdb struct')
    parser.add_argument('-p', '--pdb_dir', metavar="PDB files directory (.pdb format)", nargs='?',
                        default=sys.stdin, help='receptors directory')
    args = parser.parse_args()
    main_multi(args.pdb_dir)