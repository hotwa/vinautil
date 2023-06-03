#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :auto.py
@Description:       :自动处理脚本
@Date     :2021/09/02 15:19:15
@Author      :hotwa
@version      :1.0
https://zhuanlan.zhihu.com/p/121215784 # PyMOL 选择器的语法参考
https://blog.csdn.net/u011450367/article/details/51815130
resource: https://www.cnblogs.com/wq242424/p/13073222.html
http://blog.sina.com.cn/s/blog_7188922f0100wbz1.html
https://www.jianshu.com/p/3396e94315cb
https://blog.csdn.net/dengximo9047/article/details/101221495 # 疏水表面
http://www.mdbbs.org/thread-4064-1-1.html # 用来寻找配体口袋的残基
'''

from functools import partial
from pymol import cmd
from pathlib import Path
import pandas as pd
from json import dumps
import numpy as np
from biotools import Bdp
from loguru import logger
# from types import FunctionType,CodeType

logger.add('pymol_plugin_{time}.log')

# select ligand, resn x
# cmd.select('ligand','byres 5UEH within 6 of resn GOL resn 85P')
def autoshow(i,path,distance = 6,ionshow = False):
    """autoshow 自动展示所有配体6A以内的lines，方便查找共价键

    [extended_summary] # cmd.create('pocket',f'byres {i} within {distance} of {rawstring}') # 创建一个在当前pdb对象中配体残基周围距离为6A的口袋对象

    Arguments:
        i {[string]} -- [pdbid 例如 5EA9]

    Keyword Arguments:
        path {[string]} -- [pdb文件存放的目录，目前支持后缀为.pdb的文件，也可以在全局变量中设置好文件目录] (default: {path})
        distance {[int]} -- [显示周围原子的距离参数] (default: {6})
        ionshow {[bool]} -- [离子和有机小分子周围共价结合观察使用，默认不展示] (default: {False})

    Returns:
        [list] -- [返回ligid标识符，除去了部分离子]
    """
    cmd.reinitialize()
    p = Path(path)
    file = p.joinpath(f"{i}.pdb")
    cmd.load(file,i)
    cmd.remove('solvent metals') # 移除金属离子和溶剂
    mole = moleculeidentity(f"{i}",path)
    rawstring = 'resn ' + ' resn '.join(mole.ligIdNoion)
    print(f'{rawstring} around {distance}')
    cmd.select('ligand',f'{rawstring}')
    cmd.select('ligand_around',f'{rawstring} around {distance}') # 选择一个在当前pdb对象中配体残基周围距离为6A的口袋对象
    if ionshow: # 是否显示所有记录小分子HET条记录中的信息，对于离子与有机物显示相关共价键有效
        cmd.show('lines','ligand_part') # 显示所有HET侧链
        cmd.create('ligand_part',f'{rawstring} expand {distance}') # 单独显示小分子扩展6A周围的lines对象
    cmd.create('organ',f'organic expand {distance}')
    cmd.show('lines','organ')
    return mole.ligId

def bejson(d:dict(help='need to beauty json')) -> dict:
    return dumps(d,indent=4,ensure_ascii=False)

class covalentidentity():
    """covalentidentity 使用pymol进行识别共价键，在输出有机小分子6A以内的原子时可能出现分子的构象发生改变导致共价键距离计算有问题，还可能是与金属离子形成共价键
    """
    # __slots__ = []

    def __init__(self,pdbfilename,pdbid,path):
        """__init__ [summary]

        [extended_summary]

        Arguments:
            pdbfilename {[string]} -- pdb文件名称
            pdbid {[string]} -- 4 length id
            path {[string or pathlib obj]} -- pdb文件路径
        """
        self.pdbfilename = pdbfilename
        self.pdbid = pdbid
        self.path = Path(path)
        self.init()

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

    @classmethod
    def export_organ(cls,pdbfile,path,distance = 6,remove_water_metals = True):
        # 输出小分子结合部分的共价结合信息对象，pdb格式
        file = path.joinpath(pdbfile)
        _pdbid = covalentidentity.cleanpdb(pdbfile)
        savepath = path.joinpath(f'{_pdbid}_organ_around_{distance}_lines.pdb')
        cmd.reinitialize()
        cmd.load(filename = file,object = pdbfile)
        if remove_water_metals: cmd.remove('solvent metals') # 移除金属离子和溶剂
        cmd.create('organ',f'organic expand {distance}')
        cmd.show('lines','organ')
        cmd.save(filename=savepath,selection='organ')
        cls.connect_info_pdb_path = savepath
        return savepath

    def init(self):
        self.__create_dataframe()

    def __read_organ_pdb(self):
        exportfile = covalentidentity.export_organ(self.pdbfilename,path = Path(self.path))
        with open(file=exportfile,mode = 'r') as f:
            read_list = [i.replace('\n', '') for i in f.readlines()]
            self.read_connect = [i for i in read_list if i[:6].strip()=='CONECT']
            self.read_hetatm = [i for i in read_list if i[:6].strip()=='HETATM']
            self.read_atom = [i for i in read_list if i[:6].strip()=='ATOM']

    def __create_dataframe(self):
        self.__read_organ_pdb() # 读取初始化数据
        connect_infos = []
        for i in self.read_connect: #清洗数据
            if len(i.split()) == 3: # 策略：查找只有两个原子键的连接信息
                connect_infos.append(i)
        self.df = self.transformtodataframe(self.read_hetatm).append(self.transformtodataframe(self.read_atom))
        target_search_connect = []
        for i in connect_infos:
            l = [i for i in i.split() if i != 'CONECT']
            target_search_connect.append(l)
        self.target_search_connect = target_search_connect # 两个连接信息的列表，不一定都是共价键

    def __search_convenlent(self,distance_accuracy=1):
        true_covlent = [] # 共价键连接信息记录 有几个列表就有几个共价键链接信息
        infos = []
        locateinfos = {}
        for item in self.target_search_connect: # 对两个连接原子进行判断
            df = self.df
            # molecule chain search
            res0 = df[df['AtomNum']==item[0]]
            res1 = df[df['AtomNum']==item[1]]
            if (res0['Identifier'].all() == res1['Identifier'].all()):
                continue
            else:
                true_covlent.append(item)
                # 定位信息那个是小分子的行信息和那个是蛋白行的信息
                if res0['Identifier'].all() == 'ATOM': locateinfos['ATOM'] = res0
                if res0['Identifier'].all() == 'HETATM': locateinfos['HETATM'] = res0
                if res1['Identifier'].all() == 'ATOM': locateinfos['ATOM'] = res1
                if res1['Identifier'].all() == 'HETATM': locateinfos['HETATM'] = res1
                assert (not res0.empty and not res1.empty),(
                    'this code occur a bug'
                    'try to connet to author fix it'
                    'maybe is pymol output connect infos error'
                    'never be happen'
                )
            partinfos = self.__fill_infos(i=item,pro=locateinfos['ATOM'],mole=locateinfos['HETATM'],distance_accuracy=distance_accuracy)     
            infos.append(partinfos)
        return infos

    def convenlent_infos(self,distance_accuracy=1):
        return self.__search_convenlent(distance_accuracy=distance_accuracy)

    @property
    def cov_infos(self):
        d = self.__search_convenlent()
        return d

    @staticmethod
    def __fill_infos(i,pro,mole,distance_accuracy=1):
        distance = np.sqrt(np.square(float(pro['X'].all())-float(mole['X'].all())) + np.square(float(pro['Y'].all())-float(mole['Y'].all())) + np.square(float(pro['Z'].all())-float(mole['Z'].all())))
        infostmplate = {
            'ConnectInfos': i,
            'ResidueName':pro['ResidueName'].all(),
            'ResidueSequence':pro['ResidueSequence'].all(),
            'ChainIndentifier':pro['ChainIndentifier'].all(),
            'CovenlentAtom':{
                'protein': {
                    'AtomName':pro['AtomName'].all(),
                    'ElementSymbol':pro['ElementSymbol'].all(), 
                },
                'molecule':{
                    'AtomName':mole['AtomName'].all(),
                    'ElementSymbol':mole['ElementSymbol'].all(), 
                },
            },
            'CovenlentDistance':round(distance,distance_accuracy),
            'LigandName':mole['ResidueName'].all(),
        }
        return infostmplate

    @staticmethod
    def transformtodataframe(readlinecontent=None,first_label = 'ATOM',path = False):
        path_info = (Path(path).name.__str__(),Path(path).parent.__str__()) if path else ('unknown file name','unknown path')
        if path: readlinecontent = Bdp.read_line_list(first_column=first_label,path=path)
        if readlinecontent == None: raise ValueError('readlinecontent need')
        if first_label == ('ATOM' or 'HETATM'): # ! 注意python中的懒惰运算，及运算符特性
            # 将pdb每行读取的信息转化为dataframe 针对ATOM HETAM开头的行适用
            Identifier,AtomNum,AtomName,AtomLocationIndicator,ResidueName,ChainIndentifier,ResidueSequence,InsertionsResidue,X,Y,Z,Occupancy,Tfactor,SegmentIdentifier,ElementSymbol = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
            for i in readlinecontent:
                Identifier.append(i[:7].strip())
                AtomNum.append(i[6:11].strip())
                AtomName.append(i[12:16].strip())
                AtomLocationIndicator.append(i[16].strip())
                ResidueName.append(i[17:20].strip())
                ChainIndentifier.append(i[21].strip())
                ResidueSequence.append(i[22:26].strip())
                InsertionsResidue.append(i[26].strip())
                X.append(i[30:38].strip())
                Y.append(i[38:46].strip())
                Z.append(i[46:54].strip())
                Occupancy.append(i[54:60].strip())
                Tfactor.append(i[60:66].strip())
                SegmentIdentifier.append(i[72:76].strip())
                ElementSymbol.append(i[76:78].strip())
            df = pd.DataFrame(data={
                'Identifier':Identifier,
                'AtomNum':AtomNum,
                'AtomName':AtomName,
                'AtomLocationIndicator':AtomLocationIndicator,
                'ResidueName':ResidueName,
                'ChainIndentifier':ChainIndentifier,
                'ResidueSequence':ResidueSequence,
                'InsertionsResidue':InsertionsResidue,
                'X':X,
                'Y':Y,
                'Z':Z,
                'Occupancy':Occupancy,
                'Tfactor':Tfactor,
                'SegmentIdentifier':SegmentIdentifier,
                'ElementSymbol':ElementSymbol
            })
        elif first_label == 'LINK':
            # https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
            Record_name,Atom_name1,Alternate_location_indicator1,Residue_name1,Chain_identifier1,Residue_sequence_number1,Insertion_code1,Atom_name2,Alternate_location_indicator2,Residue_name2,Chain_identifier2,Residue_sequence_number2,Insertion_code2,Symmetry_operator_atom_1,Symmetry_operator_atom_2,Link_distance = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
            for i in readlinecontent:
                if len(i) != 78:
                    logger.exception(
                    f"[LINK label length error] not match the LINK line length\n"
                    f"[input string]:{i}\n"
                    f"check your input file:{path_info[0]} from {path_info[1]}\n"
                    "if this line not have link distance, try to caculate by covalent.bypymol method")
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
                'Record_name':Record_name,
                'Atom_name1':Atom_name1,
                'Alternate_location_indicator1':Alternate_location_indicator1,
                'Residue_name1':pd.Series(Residue_name1,dtype='object'),
                'Chain_identifier1':Chain_identifier1,
                'Residue_sequence_number1':Residue_sequence_number1,
                'Insertion_code1':Insertion_code1,
                'Atom_name2':Atom_name2,
                'Alternate_location_indicator2':Alternate_location_indicator2,
                'Residue_name2':pd.Series(Residue_name2,dtype='object'),
                'Chain_identifier2':Chain_identifier2,
                'Residue_sequence_number2':Residue_sequence_number2,
                'Insertion_code2':Insertion_code2,
                'Symmetry_operator_atom_1':Symmetry_operator_atom_1,
                'Symmetry_operator_atom_2':Symmetry_operator_atom_2,
                'Link_distance':pd.Series(Link_distance,dtype='float32'),
            })
        else:
            assert False,(
                'No return',
                'a error occurred'
            )
        return df

class color_plan():
    color1 = (('color1', '[186,182,217]'), 'purple')
    color2 = (('color2', '[233,195,153]'), 'yellow')
    color3 = (('color3', '[141,215,247]'), 'blue')
    color4 = (('color4', '[206,155,198]'), 'purple')
    color5 = (('color5', '[251,187,62]'), 'orange')
    color6 = (('color6', '[245,157,158]'), 'red')
    color7 = (('color7', '[133,188,135]'), 'green')
    colors = (color1, color2, color3, color4, color5, color6, color7)

    @staticmethod
    def defcolor():
        color_gennerate = map(lambda x:x[0],color_plan.colors)
        list(map(lambda x:cmd.set_color(*x),color_gennerate))

# error class
class PathError(BaseException):
    def __init__(self,arg):
        self.arg = arg

class moleculeidentity():
    """moleculeidentity [summary]

    [识别pdb蛋白文件中的小分子标识符]
    """

    def __init__(self,pdbfile,path):
        self.pathstr = path
        self.path = Path(path)
        self.pdbfile = pdbfile
        self._init()

    def _init(self):
        if not self.path.exists():
            raise PathError('path not exist')
        if not isinstance(self.pdbfile,str):
            raise TypeError('Pdbid must be a string!')
        if ('.pdb' not in self.pdbfile) and (len(self.pdbfile) != 4):
            raise TypeError('Pdbid must be 4 letters')
        if ('.pdb' or '.PDB') in self.pdbfile:
            raise TypeError(f'{self.pdbfile} Remove ".pdb" from input arg, add automatically')
        file_list =  list(self.path.glob('*.pdb'))
        self.path_parent = file_list[0].parent
        self.pdbfilelist = [i.name[:4].upper() for i in file_list]

    def __parse_pdb_ligid(self,ion=True):
        if self.pdbfile.upper() not in self.pdbfilelist:
            raise FileNotFoundError(f'not found {self.pdbfile} in {self.path}')
        infos_line = []
        ligId = []
        for i in self.__generate_pdb_lines():
            if self.check_line_header(i) == 'HET':
                infos_line.append(i)
        ligId = [i.split()[1] for i in infos_line]
        ligId = list(set(ligId))
        if not ion: ligId = [i for i in ligId if len(i) == 3 ] # remove ion from list
        return ligId
    
    @property
    def ligId(self):
        # return ligId include ion
        return self.__parse_pdb_ligid()

    @property
    def ligIdNoion(self):
        return self.__parse_pdb_ligid(ion=False)

    @staticmethod
    def check_line_header(line_text):
        return line_text[0:6].strip()
            
    def __generate_pdb_lines(self):
        openpdbfile = self.pdbfile + '.pdb' if '.pdb' not in self.pdbfile else self.pdbfile
        for row in open(self.path_parent.joinpath(openpdbfile),'r+'):
            yield row.strip()

class covalent(object):

    def __init__(self,pdbfilename,path,pdbid):
        self.pdbfilename = pdbfilename
        self.pdbid = covalentidentity.cleanpdb(pdbid)
        self.path_parent = Path(path)
        self.path = self.path_parent.joinpath(pdbfilename)
        self._init()

    def _init(self):
        self.link_df = covalentidentity.transformtodataframe(path = self.path,first_label='LINK')

    @classmethod
    def bypdb(cls,pdbfilename,path,pdbid):
        return cls(pdbfilename=pdbfilename,path=path,pdbid=pdbid)

    @staticmethod
    def bypymol(pdbfilename,pdbid,path):
        return covalentidentity(pdbfilename,pdbid,path)

    @property    
    def mole_id(self)->list:
        _instance=Bdp(path=self.path,sid=self.pdbid)
        return [i for i in _instance.mole_id if len(i)==3 and i != 'HOH' and i != 'SO4' and i != 'PO4'] # ! 在这里手动剔除掉硫酸根和磷酸根

    def modresdata(self)->list:
        """modresdata 标准残疾的修饰信息
        """
        modres_lst = Bdp.read_line_list(first_column='MODRES',path=self.path)
        lst = [i[12:15] for i in modres_lst] # 切片出修饰残基的小分子
        return list(set(lst))

    def get_covalent(self):
        # 从原始的pdb文件中可获取共价键信息 # ? 对其中的信息正确率为100% 来自原始的pdb文件中的信息
        res1 = self.link_df.query(f'Residue_name1 in {self.remove_modres_moleid}')
        res2 = self.link_df.query(f'Residue_name2 in {self.remove_modres_moleid}')
        df = pd.merge(res1, res2,how='outer')
        df.astype('object')
        # remove ions link with ligand
        df = df[df['Residue_name1'].str.len()==3]
        df = df[df['Residue_name2'].str.len()==3]
        return df


    @property
    def remove_modres_moleid(self)->list:
        # 清除link信息中的modres，即修饰残基的ligid
        # lst = Bdp.read_line_list(first_column='LINK',path=self.path)
        # lst = list(set(map(lambda s:s[17:20],lst)))
        clean_func = partial(self.func_dynamic_data_template,dynamic_data=self.ModifiedResidues)
        res = clean_func(origin_data=self.mole_id)
        return list(filter(lambda x:len(x.strip())==3 and x != 'HOH',res))

    @property
    def ModifiedResidues(self):
        ModifiedResiduesId = []
        for i in self.mole_id:
            if i in self.modresdata():
                ModifiedResiduesId.append(i)
        # 从SEQRES序列中获取MODRES的修饰残基的小分子，仅仅从SEQRES有时获取不全面
        seqres_lst = [i[19:].strip() for i in Bdp.read_line_list(first_column='SEQRES',path = self.path)]
        seqres_lst_string = '\r\n'.join(seqres_lst)
        seqres_lst_set = set(list(i for line in seqres_lst_string.splitlines() for i in line.split()))
        seqres_list = seqres_lst_set.intersection(set(self.mole_id))
        ModifiedResiduesId.extend(seqres_list)
        return list(set(ModifiedResiduesId))

    @staticmethod
    def func_dynamic_data_template(origin_data:list,dynamic_data:list)->list:
        # 取origin_data和dynamic_data补集
        origin_data_set = set(origin_data)
        dynamic_data_set = set(dynamic_data)
        _data = origin_data_set.difference(dynamic_data_set)
        return _data




cmd.extend('autoshow', autoshow) # pymol中自定义autoshow函数，请使用cmd.autoshow()命令

#! 下面为具体使用查找共价键的代码
# def read_link_lines(item):
#     ins1 = covalent.bypdb(pdbfilename = item,path='M:\program\\autotask\search_res1_pdb',pdbid = item[3:7])
#     res = ins1.get_covalent()
#     res['Pdb_id'] = [item[3:7] for _ in range(len(res))]
#     return res

def return_mole(series):
    # LINK 字段信息进行分析解读
    residue = [i.upper() for i in 'Gly,Ala,Val,Leu,Ile,Pro,Phe,Tyr,Trp,Ser,Thr,Cys,Met,Asn,Gln,Asp,Glu,Lys,Arg,His'.split(',')]
    if str(series[1]['Residue_name1']) in residue:
        if str(series[1]['Residue_name2']) in residue:
            # 除去两个氨基酸连接的情况
            return pd.DataFrame()
        else:
            _data =  series[1]['Residue_name2'],series[1]['Link_distance'],series[1]['Pdb_id'],series[1]['Chain_identifier2']
            _df = pd.DataFrame(data={
            'LigandName': pd.Series([_data[0]],dtype='string'),
            'Link_distance': pd.Series([_data[1]],dtype='float32'),
            'Pdb_id': pd.Series([_data[2]],dtype='string'),
            'Chain_identifier': pd.Series([_data[3]],dtype='string')
        })
            return _df
    elif str(series[1]['Residue_name2']) in residue:
        if str(series[1]['Residue_name1']) in residue:
            # 除去两个氨基酸连接的情况
            return pd.DataFrame()
        else:
            _data =  series[1]['Residue_name1'],series[1]['Link_distance'],series[1]['Pdb_id'],series[1]['Chain_identifier1']
            _df = pd.DataFrame(data={
            'LigandName': pd.Series([_data[0]],dtype='string'),
            'Link_distance': pd.Series([_data[1]],dtype='float32'),
            'Pdb_id': pd.Series([_data[2]],dtype='string'),
            'Chain_identifier': pd.Series([_data[3]],dtype='string')
        })
            return _df
    else:
        logger.exception(
            f'[MoleculeMatchError] Pdbid {series[1]["Pdb_id"]} occur a error\n'
            f'出现了两个残基名称都不属于常见的20中\n{residue}\n'
            f'source: {series[1]}'
        )
        return pd.DataFrame() # 两个小分子连接的情况，去除

if __name__ == '__main__':
    ...