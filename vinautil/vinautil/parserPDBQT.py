#!/usr/bin/ehon
# -*- encoding: utf-8 -*-
'''
@file        :parserPDBQT.py
@Description:       : parser PDBQT file from autodock vina
@Date     :2022/10/07 17:08:03
@Author      :hotwa
@version      :1.0

add 2022 10 07 SumNoHAtom
'''
# from typing import TypeVar, Text
import re
from pathlib import Path
# from io import StringIO,BytesIO
from pdbparse import Bdp

class qtparser(object):
    __slots__ = ["file", "path", "file_stem", "file_obj", "lines_num"]

    def __init__(self,file):
        self.file = file
        self.path = Path(file).parent
        self.file_stem = Path(file).stem
        self.file_obj = Path(file)
        self.init()

    def init(self):
        with open(self.file, 'r') as f:
            self.lines_num = len(f.readlines())

    @property
    def module_number(self):
        return len(self.getmodules())

    def getmodules(self):
        module_lst = []
        module = []
        f = open(self.file, 'r')
        for i in range(self.lines_num):
            line = f.readline()
            if line.strip() != 'ENDMDL':
                module.append(line)
            else:
                module.append(line)
                module_lst.append(module)
                module = []
        return module_lst

    def extract_pose(self, sequence:int, out_file:Path=None, overwrite:bool=True):
        """
        extract pdbqt pose
        :param file: output file
        :param sequence: extract pose number
        :return:
        """
        lst = self.getmodules()
        _string = ''.join(lst[sequence])
        if out_file:
            _p = Path(out_file)
            if not _p.parent.exists(): _p.parent.mkdir(parents=True)
            with open(out_file, 'w', encoding='utf-8') as f:
                f.write(_string)
            return out_file
        else:
            return _string

    @classmethod
    def write_pose(cls, file):
        self = cls(file)
        lst = self.getmodules()
        for i in range(len(lst)):
            _string = ''.join(lst[i])
            _p = self.path.joinpath(f'{self.file_stem}')
            if not _p.exists(): _p.mkdir(parents = True)
            with open(_p.joinpath(self.file_stem + f'_{i}.pdbqt').__str__(), 'w') as f:
                f.write(_string)
        return self

    def scores(self):# get affinity from pdbqt file
        score_lst = []
        module_lst = self.getmodules()
        for i,j in enumerate(module_lst):
            # print(f'MODEL {i+1} affinity')
            for ii in j:
                if 'VINA RESULT' in ii:
                    # print(f'locate line {ii}')
                    score_lst.append(ii.split()[-3])
        return score_lst

    @property
    def energies(self):
        return self.scores()
    
    # 通过pdbqt文件计算非H原子格式
    def read_file(self):
        """read_file read pdbqt file as a dict, key is "MODEL N", value is lines.

        [extended_summary]

        Arguments:
            file {[type]} -- [description]

        Returns:
            [type] -- [description]
        """
        d,s,model = {},'',None
        f = open(self.file, 'r')
        while contentline := f.readline():
            if 'MODEL' in contentline:
                model = contentline.strip()
            else:
                j = bool(contentline.strip() == 'ENDMDL')
                if j and model:
                    d[model] = s
                    s = ''
                if not j: s += contentline
        f.close()
        return d

    def SumNoHAtom(self, structure_id = 'MODEL 1'): # 计算非H原子个数
        d = self.read_file()
        model = d[structure_id].split('\n')
        # clean data
        model = [i for i in model if i]
        model = [i for i in model if i.split()[0] == 'HETATM']
        df = Bdp.transformtodataframe(readlinecontent = model, first_label = 'HETATM', keep_space=False)
        atom_lst = [re.match(string = i, pattern='\D*').group(0) for i in df['AtomName'].to_list()]
        n = 0
        for i in atom_lst:
            if 'H' == i:
                ...
            else:
                n += 1
        return n

