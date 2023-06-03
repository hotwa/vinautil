#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :pdbparse.py
@Description:       : hanlde PDB file
@Date     :2022/10/13 15:28:49 update
@Author      :hotwa
@version      :1.0
'''
# need biopython
import os, time
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO, Select, parse_pdb_header, PDBList
from pathlib import Path

try:
    from .typecheck import typeassert
except (ImportError, ModuleNotFoundError, AttributeError, NameError, SyntaxError):
    from typecheck import typeassert

file_dir = os.path.dirname(os.path.realpath(__file__))  # 当前python文件的绝对路径
path_obj = Path(file_dir)


# timestring = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# local import data
def get_bonding(file='amino_bonding_atom.csv'):
    file = path_obj.joinpath(file)
    return pd.read_csv(file)


def get_complex(file='./new_complex_finish_update_pdb_site.xlsx', same_chain=False):
    """
    :param file: xlsx file
    :param same_chain: Chain_identifier_mole == Residue_chain
    :return:
    """
    file = path_obj.joinpath(file)
    df = pd.read_excel(file)
    if same_chain:
        return df[df['Chain_identifier_mole'] == df['Residue_chain']]
    else:
        return df


class PDBLengthError(BaseException):
    def __init__(self, arg):
        self.arg = arg


class ChainSelect(Select):
    def __init__(self, chain_string='A'):
        super(Select, self).__init__()
        self.chain_string = chain_string

    def accept_chain(self, chain):
        if str(chain.__repr__()) == '<Chain id={}>'.format(self.chain_string):  # judge chain name
            return 1
        else:
            return 0


class ResidueSelect(Select):
    """ResidueSelect ues by select

    [extended_summary]

    :param Select: [description]
    :type Select: [type]
    """

    def __init__(self, residue):
        super(Select, self).__init__()
        self.residue = residue

    def accept_residue(self, residue):
        if residue.get_resname() == self.residue:
            return 1
        else:
            return 0


# ! error class for Bdp
class DoiError(BaseException):
    def __init__(self, arg):
        self.arg = arg


class ChainError(BaseException):
    def __init__(self, arg):
        self.arg = arg


@typeassert(path=Path, sid=str)
class Bdp(object):
    __slots__ = ['path', 'sid', 'mole_id', 'mole_struct', 'create_path', '__dict__']

    def __init__(self, path: Path, sid: str):
        self.path = path
        self.create_path = path.parent
        self.sid = sid
        self.mole_id = []
        self.mole_struct = []
        # 初始化数据
        self.read_formula()

    @property
    def chainnum(self):
        s = self.read_struct()
        first_model = s[0]
        return len(first_model)

    @property
    def chainstr(self):
        return self.getchain_str()

    @property
    def headerinfos(self):
        return parse_pdb_header(self.path)

    @property
    def doi(self):
        _journal = self.headerinfos['journal']
        _doi = _journal.split()[-1]
        if _doi.startswith('10.'):
            return _doi
        else:
            raise DoiError(f'current pdb does not doi, {_journal}')

    def residuesequence(self, label):  # 'HETATM' 'ATOM'
        df = self.transformtodataframe(path=self.path, first_label=label)
        return set(df['ResidueSequence'].to_list())

    def value_check(self):
        try:
            if len(self.mole_id) == len(self.mole_struct):
                print('check pass!')
            else:
                raise ValueError('molecule identifier do not match formula identifier! manual check')
        except:
            raise ValueError

    @staticmethod
    def read_line_list(first_column, path):
        """read_line_list 读取第一列为first_column字符串的行，存为列表返回

        [extended_summary]

        Arguments:
            first_column {[string]} -- [pdb文件内容第一列常常为大写]
            path {[string]} -- [路径]

        Returns:
            [iter] -- [含有first_column字符串的生成器]
        """
        stringiter = Bdp.stringlinesiter(file=path)
        stringiter = map(lambda x: x.strip(), stringiter)
        stringiter = filter(lambda x: bool(x), stringiter)  # clean empty
        return filter(lambda x: x.split()[0] == first_column, stringiter)

    @staticmethod
    def stringlinesiter(file):
        with open(file, 'r+', encoding='utf-8') as f:
            yield from f.readlines()

    def mole_site(self, lig: str, chain: str) -> int:
        readmolelist = self.split_molecule(residue=lig, chain=chain)
        with open(readmolelist, 'r'):
            df = Bdp.transformtodataframe(first_label='HETATM', path=readmolelist)
        return df['ResidueSequence'].tolist()[0]

    def _split_molecule_from_chainlist(self, chainlist: list, s: object, residue: str) -> list:
        # 对每条链上都与这个小分子结合的链上的结合信息都提取出来
        save_list, remove_list = [], []
        for i in chainlist:
            io = PDBIO()
            io.set_structure(s[0][i])
            sidc = self.sid.replace('.pdb', '') if '.pdb' in self.sid else self.sid
            savename = f'{sidc}_{i}_{residue}.pdb'
            create_dir = self.create_path.joinpath('split_molecule')
            if not create_dir.exists(): create_dir.mkdir()
            path = create_dir.joinpath(savename)
            save_list.append(path.__str__())
            io.save(path.__str__(), ResidueSelect(residue))
            if path.exists():
                with open(path, 'r+') as f:
                    content = f.readline()
                    if content.strip() == 'END':
                        # print(f'{savename} this chain have not molecule remove it')
                        remove_list.append(path)  # 该链没有小分子，删除文件
        for i in remove_list:
            os.remove(i)  # remove the empty molecule file
        header_info = f'''---
title: Split molecule from PDB
date: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}
type: pdb-{self.sid}
---
[TOC]
# Description
split {','.join(self.mole_id)}
'''
        create_path_md = self.create_path.joinpath('split_molecule', 'README.md')
        with open(create_path_md, 'w+') as f:
            f.write(header_info)
        return save_list

    def split_molecule(self, residue: str, chain: str = None,
                       keep_one_struct: bool = True) -> list:  # ! 对于修饰残基的提取不完美 bug
        '''
        split_molecule 从pdb文件中提取小分子
        :param residue: 小分子的名字，ligandID
        :param chain: 小分子所在的链，如果为None则提取所有链
        :return: 返回小分子的路径 list
        # ! 需要对提取的小分子构象二次检查，小分子构象可能有两个 residue_keep_struct
        '''
        residue = str(residue)
        s = self.read_struct()
        chainlist = self.getchain_str()  # 该蛋白获取所有链id编号
        if chain == None:
            # 提取所有链上的小分子
            save_list = self._split_molecule_from_chainlist(chainlist=chainlist, s=s, residue=residue)
            save_list = list(map(Bdp.residue_keep_struct_file, save_list))
        else:
            if chain in chainlist:
                chainlist = [chain]
            else:
                raise ValueError(f'chain ID :{chain} not in {self.sid}, file: {self.path}')
            save_list = self._split_molecule_from_chainlist(chainlist, s, residue)
        save_list = list(map(Bdp.residue_keep_struct_file, save_list))
        return save_list

    @staticmethod
    def residue_keep_struct_file(file):  # design for split_molecule method
        # smart to recognize the chain
        df_mol = Bdp.transformtodataframe(first_label='HETATM', path=file)
        try:
            construct_id = Bdp.get_residue_AtomLocationIndicator(df_mol)
            id = construct_id[0]
        except IndexError as e:
            raise e(f'muti construct get failed, please check {file}')
        df_item = Bdp.residue_keep_struct(file, id=id, first_label='HETATM')
        if df_item.empty:
            raise ValueError(f'error in {file} to get residue_keep_struct')
        else:
            out_file = Bdp.dataframe2pdb(df_item, file)
            if Path(out_file).exists() and out_file == file:
                return out_file  # transform success

    @staticmethod
    def clean_struct(file, label):  # 清除所有多重构象的结构构象 AtomLocationIndicator 中保留第一个
        df0 = Bdp.transformtodataframe(path=file, first_label=label, keep_space=False)
        seq_lst = list(set(df0['ResidueSequence'].to_list()))
        df1 = Bdp.transformtodataframe(first_label=label, path=file, keep_space=True)
        df = Bdp.transformtodataframe(first_label=label, path=file, keep_space=False)
        # 搜集含有多重构象的残基序列,并且保留第一个构象
        remove_index_lst = []  # 需要删除的行(索引)
        for site in seq_lst:
            df_muta = df[df['ResidueSequence'] == str(site)]  # 定位残基
            construct_id = Bdp.get_residue_AtomLocationIndicator(df_muta)
            if construct_id:  # true 有多重构象
                id = construct_id[0]
                # 删除除了id构象的其他构象信息
                all_index = df_muta.index
                keep_index = df_muta[df_muta['AtomLocationIndicator'] == id].index  # 需要保留的侧链信息
                public_index = df_muta[df_muta['AtomLocationIndicator'] == ''].index  # 主链信息
                remove_index = list(set(all_index) - set(keep_index) - set(public_index))  # 需要删除的索引信息
                remove_index_lst.extend(remove_index)
        remove_index_data = pd.Index(data=remove_index_lst
                                     , dtype='int64')
        df_out = df1.drop(index=remove_index_data, axis=0)
        df_out.index = [i for i in range(len(df_out))]
        return df_out  # 使用 Bdp.dataframe2pdb 转换为pdb文件

    @staticmethod
    def keep_residue(file, site, first_label: str = 'ATOM', id='auto'):
        """
        keep_residue 智能保留残基
        :param file: read file, pdb file
        :param site: residue site
        :param first_label: 'ATOM' or 'HETATM'
        :param id: keep construct id(like A, B……) default: first id(A)
        :return: pd.DataFrame
        """
        df1 = Bdp.transformtodataframe(first_label, path=file, keep_space=True)
        df = Bdp.transformtodataframe(first_label, path=file, keep_space=False)
        # 清除原始蛋白中的site位点信息
        df_remove = df[df['ResidueSequence'] != str(site)]
        # 定位残基
        df_muta = df[df['ResidueSequence'] == str(site)]
        df_muta_num = len(df_muta)
        construct_id = Bdp.get_residue_AtomLocationIndicator(df_muta)
        if construct_id:
            if id == 'auto':
                id = construct_id[0]
            else:
                if id not in construct_id:
                    raise ValueError(f'{file} 在 {site} 没有 {id} 构象信息, 可选构象为 {construct_id}')
                id = id.upper()
        else:
            raise ValueError(f'{file} 在 {site} 没有构象信息')
        # 删除除了id构象的其他构象信息
        all_index = df_muta.index
        keep_index = df_muta[df_muta['AtomLocationIndicator'] == id].index  # 需要保留的侧链信息
        public_index = df_muta[df_muta['AtomLocationIndicator'] == ''].index  # 主链信息
        remove_index = list(set(all_index) - set(keep_index) - set(public_index))
        remove_index = pd.Index(data=remove_index, dtype='int64')
        df_out = df1.drop(index=remove_index, axis=0)
        df_out.index = [i for i in range(len(df_out))]
        return df_out

    @staticmethod
    def get_residue_AtomLocationIndicator(df):
        # 返回残基的DataFrame构象的标识符号，没有则返回空列表
        lst = [i for i in list(set(df['AtomLocationIndicator'].tolist())) if i]
        lst.sort()
        return lst

    @staticmethod
    def residue_keep_struct(file, id: str = 'A', first_label: str = 'ATOM'):
        # 需要重新修改该函数，针对摸个特定的残基进行保留相应的构象，而不是同一保留所有的构象，例如A构象
        # 有些残基在A构象中不存在，但是在B构象中存在，会导致报错，残基会出现缺失的情况。
        '''
        residue_keep_construct 保留构象，删除其他构象
        :param file: read file, pdb file
        :param id: keep construct id（like A，B……）
        :param first_label: ‘’ATOM‘’ or ‘’HETATM‘’
        :return: pd.DataFrame
        '''
        df = Bdp.transformtodataframe(first_label, path=file, keep_space=True)
        remove_index = []
        for i in df.iterrows():
            atomsign = i[1]['AtomLocationIndicator']
            if bool(atomsign.strip()):
                if atomsign == id.upper():
                    df.loc[i[0], 'AtomLocationIndicator'] = ' '
                else:
                    print(atomsign)
                    remove_index.append(i[0])
        df = df.drop(index=df.index[remove_index], axis=0)
        df.index = [i for i in range(len(df))]
        return df

    def getchain_str(self):
        """getchain_str retrun chain symbol name

        [extended_summary]

        :return: [description]
        :rtype: list
        """
        _l = []
        p = PDBParser(PERMISSIVE=1)
        s = p.get_structure(self.sid, self.path)
        chain_gennerate = s[0].copy().get_chains()
        for i in chain_gennerate:
            _l.append(i.id)
            _l.sort()
        return _l

    def site_exists(self, residue, num):
        """site_exists 判断改蛋白的某个位点是否存在，并返回该位点所在的链

        用于判断结合位点的链在哪里

        :param residue: 残基名称
        :type residue: string
        :param num: 残基编号
        :type num: int
        :return: chain
        :rtype: list
        """
        _l = []
        num = int(num)
        residue = residue.strip().upper()
        chainlist = self.getchain_str()
        s = self.read_struct()
        for i in chainlist:
            try:
                resname_from_pdb = s[0][i][num].get_resname()
                if residue == resname_from_pdb:
                    _l.append(i)
            except KeyError as e:
                print(e, f'chain {i} not exist {residue}-{num}')
        return _l

    def link_info(self, qurey_key=None):
        dfl = self.transformtodataframe(first_label='LINK', path=self.path)
        if qurey_key == None:
            return dfl
        else:
            return dfl[(dfl['Residue_name1'] == qurey_key) | (dfl['Residue_name2'] == qurey_key)]

    def link_mole_site(self, mole_name: str):
        mole_name = mole_name.upper()
        df = self.link_info(qurey_key=mole_name)
        if df.empty:
            return None
        elif len(df) == 1:
            item = df.values.tolist()[0]
            if item[3] == mole_name:
                return item[5]
            elif item[9] == mole_name:
                return item[11]
            else:
                raise ValueError('impossible')
        else:
            sites_lst = []
            data = df.values.tolist()
            for item in data:
                if item[3] == mole_name:
                    sites_lst.append(item[5])
                elif item[9] == mole_name:
                    sites_lst.append(item[11])
                else:
                    ...
            assert len(sites_lst) == len(data)
            return sites_lst

    @classmethod
    def split_chains(cls, file: Path, sid, out_filename=None, extract='auto', clean=True):
        directory = Path(file).parent
        if not directory.exists(): directory.mkdir(parents=True)
        self = cls(file, sid)
        file = cls.clean_pdb(file) if clean else file  # keep ATOM line
        auth_chainlist = self.getchain_str()
        p = PDBParser(PERMISSIVE=1)
        s = p.get_structure(self.sid, file)
        file = Path(file)
        if extract == 'auto':
            for j in chainlist:
                filename = f'{file.stem}_{j}.pdb'
                self._split_chain_save_struct(obj=s[0][j], sid=self.sid, file=filename)
        elif isinstance(extract, list):
            for j in extract:
                if j in auth_chainlist:
                    filename = f'{file.stem}_{j}.pdb'
                    self._split_chain_save_struct(obj=s[0][j], sid=self.sid, file=filename)
                else:
                    print(f'chain {j} not exist')
        elif isinstance(extract, str):
            if extract in auth_chainlist:
                filename = out_filename if out_filename else f'{file.stem}_{extract}.pdb'
                self._split_chain_save_struct(obj=s[0][extract], sid=self.sid, file=filename)
            else:
                print(f'chain {extract} not exist')
        else:
            ...
        return self

    def _split_chain_save_struct(self, obj, sid, file):
        io = PDBIO()
        io.set_structure(obj)
        sidc = sid.replace('.pdb', '') if '.pdb' in sid else sid
        if not Path(file).parent.exists(): Path(file).parent.mkdir(parents=True)
        io.save(file.__str__(), ChainSelect(obj.id))
        return file

    @staticmethod
    def clean_pdb(file, outfile=None):
        """clean_pdb 清除杂原子，保留ATOM

        _extended_summary_

        Keyword Arguments:
            outfile {string} -- _description_ (default: {None})
            path {path} -- _description_ (default: {None})

        Returns:
            object -- path
        """
        __lst = Bdp.read_line_list('ATOM', file)
        _f = Path(outfile) if outfile else Path(file).parent.joinpath(f"{Path(file).stem}.clean.pdb")
        if _f.exists(): _f.unlink()
        with open(_f.__str__(), 'w+') as f:
            for i in __lst:
                f.write(f"{i}\n")
        return _f

    def read_struct(self):
        """read_struct read pdb structure

        https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr11.html#id3
        一个 Structure 对象的整体布局遵循称为SMCRA（Structure/Model/Chain/Residue/Atom，结构/模型/链/残基/原子）的体系架构：

        结构由模型组成
        模型由多条链组成
        链由残基组成
        多个原子构成残基
        Returns:
            [obj] -- [Model]
        """
        p = PDBParser(PERMISSIVE=1)
        s = p.get_structure(self.sid, self.path)
        return s

    def get_chain(self, chain_str):
        s = self.read_struct()
        return s[0][chain_str]

    def read_formula(self):
        """read_formula read pdb file molecule information and identifier in this pdb file(just like residue)

        [extended_summary]

        :raises OSError: Do not find this file
        """
        _path = self.path
        try:
            with open(_path, 'r') as f:
                text_line = f.readlines()
                _l = []
                for i in text_line:
                    li = i.split()
                    if li:
                        if li[0] == 'FORMUL':
                            _l.append(li)
                            self.mole_id.append(li[2])
                for ii in _l:
                    self.mole_struct.append(''.join(ii[3:]))
        except Exception as e:
            raise e

    def mole_link():
        Bdp.transformtodataframe(first_label='LINK', )

    @staticmethod
    def transformtodataframe(first_label='ATOM', readlinecontent=None, path=False, keep_space: bool = True):
        # https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
        path_info = (Path(path).name.__str__(), Path(path).parent.__str__()) if path else (
            'unknown file name', 'unknown path')
        if path: readlinecontent = Bdp.read_line_list(first_column=first_label, path=path)
        if readlinecontent == None: raise ValueError('readlinecontent need')
        if first_label == 'ATOM' or first_label == 'HETATM':  # ! 注意python中的懒惰运算，及运算符特性
            # 将pdb每行读取的信息转化为dataframe 针对ATOM HETAM开头的行适用
            Identifier, AtomNum, AtomName, AtomLocationIndicator, ResidueName, ChainIndentifier, ResidueSequence, InsertionsResidue, X, Y, Z, Occupancy, Tfactor, SegmentIdentifier, ElementSymbol = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
            for i in readlinecontent:
                if keep_space:  # for restore to pdb file
                    Identifier.append(i[:7])
                    AtomNum.append(i[6:12])  # 这里应该是切片6:12，但是PDB官网是6:11，导致AtomNum的值可能是错的，这里预留了一位数字，在超过10个原子时会出现两位数的情况
                    AtomName.append(i[12:16])
                    AtomLocationIndicator.append(i[16])
                    ResidueName.append(i[17:20])
                    ChainIndentifier.append(i[21])
                    ResidueSequence.append(i[22:26])
                    InsertionsResidue.append(i[26])
                    X.append(i[30:38])
                    Y.append(i[38:46])
                    Z.append(i[46:54])
                    Occupancy.append(i[54:60])
                    Tfactor.append(i[60:66])
                    SegmentIdentifier.append(i[72:76])
                    ElementSymbol.append(i[76:78])
                else:
                    Identifier.append(i[:7].strip())
                    AtomNum.append(
                        i[6:12].strip())  # 这里应该是切片6:12，但是PDB官网是6:11，导致AtomNum的值可能是错的，这里预留了一位数字，在超过10个原子时会出现两位数的情况
                    AtomName.append(i[12:16].strip())
                    AtomLocationIndicator.append(i[16].strip())
                    ResidueName.append(i[17:20].strip())
                    ChainIndentifier.append(i[21].strip())
                    ResidueSequence.append(i[22:26].strip())
                    InsertionsResidue.append(i[26].strip())
                    X.append(float(i[30:38].strip()))
                    Y.append(float(i[38:46].strip()))
                    Z.append(float(i[46:54].strip()))
                    Occupancy.append(i[54:60].strip())
                    Tfactor.append(i[60:66].strip())
                    SegmentIdentifier.append(i[72:76].strip())
                    ElementSymbol.append(i[76:78].strip())
            df = pd.DataFrame(data={
                'Identifier': Identifier,
                'AtomNum': AtomNum,
                'AtomName': AtomName,
                'AtomLocationIndicator': AtomLocationIndicator,
                'ResidueName': ResidueName,
                'ChainIndentifier': ChainIndentifier,
                'ResidueSequence': ResidueSequence,
                'InsertionsResidue': InsertionsResidue,
                'X': X,
                'Y': Y,
                'Z': Z,
                'Occupancy': Occupancy,
                'Tfactor': Tfactor,
                'SegmentIdentifier': SegmentIdentifier,
                'ElementSymbol': ElementSymbol
            })
        elif first_label == 'LINK':
            # https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
            Record_name, Atom_name1, Alternate_location_indicator1, Residue_name1, Chain_identifier1, Residue_sequence_number1, Insertion_code1, Atom_name2, Alternate_location_indicator2, Residue_name2, Chain_identifier2, Residue_sequence_number2, Insertion_code2, Symmetry_operator_atom_1, Symmetry_operator_atom_2, Link_distance = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
            for i in readlinecontent:
                if len(i.strip()) != 78:
                    # logger.exception(
                    # f"[LINK label length error] not match the LINK line length\n"
                    # f"[input string]:{i}\n"
                    # f"check your input file:{path_info[0]} from {path_info[1]}\n"
                    # "if this line not have link distance, try to caculate by covalent.bypymol method")
                    pass
                else:
                    Record_name.append(i[:6].strip())
                    Atom_name1.append(i[12:16].strip())
                    Alternate_location_indicator1.append(i[16].strip())
                    Residue_name1.append(i[17:20].strip())
                    Chain_identifier1.append(i[21].strip())
                    Residue_sequence_number1.append(i[22:26].strip())
                    Insertion_code1.append(i[26].strip())
                    Atom_name2.append(i[42:46].strip())
                    Alternate_location_indicator2.append(i[46].strip())
                    Residue_name2.append(i[47:50].strip())
                    Chain_identifier2.append(i[51].strip())
                    Residue_sequence_number2.append(i[52:56].strip())
                    Insertion_code2.append(i[56].strip())
                    Symmetry_operator_atom_1.append(i[59:65].strip())
                    Symmetry_operator_atom_2.append(i[66:72].strip())
                    Link_distance.append(i[73:78].strip())
            df = pd.DataFrame(data={
                'Record_name': Record_name,
                'Atom_name1': Atom_name1,
                'Alternate_location_indicator1': Alternate_location_indicator1,
                'Residue_name1': pd.Series(Residue_name1, dtype='object'),
                'Chain_identifier1': Chain_identifier1,
                'Residue_sequence_number1': Residue_sequence_number1,
                'Insertion_code1': Insertion_code1,
                'Atom_name2': Atom_name2,
                'Alternate_location_indicator2': Alternate_location_indicator2,
                'Residue_name2': pd.Series(Residue_name2, dtype='object'),
                'Chain_identifier2': Chain_identifier2,
                'Residue_sequence_number2': Residue_sequence_number2,
                'Insertion_code2': Insertion_code2,
                'Symmetry_operator_atom_1': Symmetry_operator_atom_1,
                'Symmetry_operator_atom_2': Symmetry_operator_atom_2,
                'Link_distance': pd.Series(Link_distance, dtype='float32'),
            })
        else:
            assert False, (
                'No return',
                'a error occurred'
            )
        return df

    @staticmethod
    def dataframe2pdb(df, out_file):
        """
        将dataframe转化为pdb文件, support ATOM type: ATOM, HETATM
        """
        with open(out_file, 'w') as f:
            empyt_line = list(' ' * 80)
            for i in df.iterrows():
                if i[1]['Identifier'].strip() == 'ATOM' or i[1]['Identifier'].strip() == 'HETATM':
                    c = i[1].tolist()
                    empyt_line[:7] = list(c[0])
                    empyt_line[6:12] = list(c[1])
                    empyt_line[12:16] = list(c[2])
                    empyt_line[16] = list(c[3])
                    empyt_line[17:20] = list(c[4])
                    empyt_line[21] = list(c[5])
                    empyt_line[22:26] = list(c[6])
                    empyt_line[26] = list(c[7])
                    empyt_line[30:38] = list(c[8])
                    empyt_line[38:46] = list(c[9])
                    empyt_line[46:54] = list(c[10])
                    empyt_line[54:60] = list(c[11])
                    empyt_line[60:66] = list(c[12])
                    empyt_line[72:76] = list(c[13])
                    empyt_line[76:78] = list(c[14])
                    empyt_line[78:80] = [' ', ' ']
                    for i, j in enumerate(empyt_line):
                        if isinstance(j, list) and len(j) == 1:
                            empyt_line[i] = j[0]
                    cs = ''.join(empyt_line)
                    f.writelines(cs + '\n')
                    empyt_line = list(' ' * 80)
        return out_file

    @staticmethod
    def cleanpdb(before_string):
        """
        处理传入参数pdb 例如: pdb1b12.ent 1b12.pdb 转化成 1b12
        """
        # print(f'try to format {before_string}')
        if '.ent' in before_string:
            before_string = before_string[3:-4].lower()
        elif '.pdb' in before_string:
            before_string = before_string[:4].lower()
        elif len(before_string) == 4:
            before_string = before_string.lower()
        else:
            if len(before_string) != 4:
                raise ValueError(f'length out of 4 {before_string}')
        return before_string

    def modresdata(self) -> list:
        """modresdata 标准残疾的修饰信息
        # 从SEQRES序列中获取MODRES的修饰残基的小分子，仅仅从SEQRES有时获取不全面
        使用ModifiedResidues更加全面
        """
        modres_lst = self.read_line_list(first_column='MODRES', path=self.path)
        lst = [i[12:15] for i in modres_lst]  # 切片出修饰残基的小分子
        return list(set(lst))

    def get_covalent(self):
        # 从原始的pdb文件中可获取共价键信息 # ? 对其中的信息正确率为100% 来自原始的pdb文件中的信息
        res1 = self.link_df.query(f'Residue_name1 in {self.remove_modres_moleid}')
        res2 = self.link_df.query(f'Residue_name2 in {self.remove_modres_moleid}')
        df = pd.merge(res1, res2, how='outer')
        df.astype('object')
        # remove ions link with ligand # ! some molecule like 5H
        # df = df[df['Residue_name1'].str.len() == 3]
        # df = df[df['Residue_name2'].str.len() == 3]
        return df

    @property
    def ModifiedResidues(self):
        ModifiedResiduesId = []
        for i in self.mole_id:
            if i in self.modresdata():
                ModifiedResiduesId.append(i)
        # 从SEQRES序列中获取MODRES的修饰残基的小分子，仅仅从SEQRES有时获取不全面
        seqres_lst = [i[19:].strip() for i in Bdp.read_line_list(first_column='SEQRES', path=self.path)]
        seqres_lst_string = '\r\n'.join(seqres_lst)
        seqres_lst_set = set(list(i for line in seqres_lst_string.splitlines() for i in line.split()))
        seqres_list = seqres_lst_set.intersection(set(self.mole_id))
        ModifiedResiduesId.extend(seqres_list)
        return list(set(ModifiedResiduesId))

    @property
    def remove_modres_moleid(self) -> list:
        # 清除link信息中的modres，即修饰残基的ligid
        # lst = Bdp.read_line_list(first_column='LINK',path=self.path)
        # lst = list(set(map(lambda s:s[17:20],lst)))
        clean_func = partial(self.func_dynamic_data_template, dynamic_data=self.ModifiedResidues)
        res = clean_func(origin_data=self.mole_id)
        return list(filter(lambda x: len(x.strip()) == 3 and x != 'HOH', res))

    @staticmethod
    def func_dynamic_data_template(origin_data: list, dynamic_data: list) -> list:
        # 取origin_data和dynamic_data补集
        origin_data_set = set(origin_data)
        dynamic_data_set = set(dynamic_data)
        _data = origin_data_set.difference(dynamic_data_set)
        return _data