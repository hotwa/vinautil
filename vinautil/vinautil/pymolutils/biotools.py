from email import header
from unittest.mock import Base


#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :biotools.py
@Description:       :
@Date     :2021/06/22 16:36:34
@Author      :hotwa
@version      :1.0
'''
import os,re,time
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO,Select,parse_pdb_header,PDBList
from tools import mdcheck,get_path_contents
import requests
import random,string
from pathlib import Path

file_dir = os.path.dirname(os.path.realpath(__file__)) # 当前python文件的绝对路径
residue = [i.upper() for i in 'Gly,Ala,Val,Leu,Ile,Pro,Phe,Tyr,Trp,Ser,Thr,Cys,Met,Asn,Gln,Asp,Glu,Lys,Arg,His'.split(',')]

""" getpdbfile function"""
agent = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.149 Safari/537.36 -- '+''.join(random.sample(string.ascii_letters+string.digits, 32))
HEADER = {'User-Agent': agent}

def getpdbfile(pdbid,dir): # 下载备用,不是很稳定
    pdbid = pdbid.upper()
    res = os.path.exists(f'{pdbid}.pdb')
    if not res:
        r = requests.get(url = f'https://files.rcsb.org/download/{pdbid}.pdb1',verify=False,headers=HEADER)
        _file_path = os.path.join(dir,f'{pdbid}.pdb')
        with open(_file_path,'wb+') as f:
            f.write(r.content)

""" getpdbfile function"""


class PDBLengthError(BaseException):
    def __init__(self,arg):
        self.arg = arg

def downloadpdbfile(pdbid,pdir,overwrite=True,changename=True):
    """downloadpdbfile [summary]
    
    [extended_summary]
    
    :param pdbid: like code 5xv7
    :type pdbid: string
    :param pdir: save pdb file path
    :type pdir: string
    :param overwrite: overwrite pdb file, defaults to True
    :type overwrite: bool, optional
    :param changename: change .ent to .pdb, defaults to True
    :type changename: bool, optional
    """
    if isinstance(pdbid,list):
        for i in pdbid:
            downloadpdbfile_sub(pdbid=i,dir = pdir,overwrite=overwrite,changename=changename)
    elif isinstance(pdbid,str):
        if len(pdbid) == 4:
            downloadpdbfile_sub(pdbid=pdbid,dir = pdir,overwrite=overwrite,changename=changename)
    else:
        raise PDBLengthError('pdbid length error!')


def downloadpdbfile_sub(pdbid,dir,overwrite,changename):
    _pdbname = f'pdb{pdbid}.ent'
    _pdbname_new = f'{pdbid}.pdb'
    if not overwrite: 
        _path = os.path.join(dir,_pdbname_new)
        if os.path.exists(_path): 
            ...
        else:
            pdbl = PDBList()
            try:
                pdbl.retrieve_pdb_file(pdbid,file_format='pdb',pdir = dir,overwrite=overwrite)
            except FileNotFoundError: # 调用biopython下载接口失败，使用自己的接口
                getpdbfile(pdbid,dir)
            if changename:
                srcFile = os.path.join(dir,_pdbname)
                dstFile = os.path.join(dir,_pdbname_new)
                try:
                    os.rename(srcFile,dstFile)
                except FileExistsError as e: # 如果重名名的文件存在，则表示该文件已经下载
                    ...
                except Exception as e:
                    raise e
    elif overwrite: # 覆盖之前已经下载的pdb文件
        pdbl = PDBList()
        try:
            pdbl.retrieve_pdb_file(pdbid,file_format='pdb',pdir = dir,overwrite=overwrite)
        except FileNotFoundError: # 调用biopython下载接口失败，使用自己的接口
            getpdbfile(pdbid,dir)
        if changename:
            srcFile = os.path.join(dir,_pdbname)
            dstFile = os.path.join(dir,_pdbname_new)
            try:
                os.rename(srcFile,dstFile)
            except Exception as e:
                raise e
    else:
        print('error overwrite params!')

class ChainSelect(Select):
    def __init__(self,chain_string='A'):
        super(Select,self).__init__()
        self.chain_string = chain_string

    def accept_chain(self, chain):
        if str(chain.__repr__()) == '<Chain id={}>'.format(self.chain_string): # judge chain name
            return 1
        else:
            return 0

class ResidueSelect(Select):
    """ResidueSelect ues by select 
    
    [extended_summary]
    
    :param Select: [description]
    :type Select: [type]
    """
    def __init__(self,residue):
        super(Select,self).__init__()
        self.residue = residue

    def accept_residue(self, residue):
        if residue.get_resname()==self.residue:
            return 1
        else:
            return 0

def return_path(path):
    # 判断是否为绝对路径
    path_split_string = os.path.split(path)
    if os.path.isdir(path_split_string[0]) and os.path.isfile(path):
        return path # 输入绝对路径，返回绝对路径
    elif bool(path_split_string[0]) == False:
        return os.path.join(file_dir,path) # 尝试在当前目录下自动补全路径
    elif  not (path_split_string[1][-4:] == '.pdb' or '.ent'):
        raise ValueError(f'this file, {path} is not a valid pdb file')
    else:
        raise ValueError(f'valid path: {path}')

# ! error class for Bdp
class DoiError(BaseException):
    def __init__(self,arg):
        self.arg = arg

class ChainError(BaseException):
    def __init__(self,arg):
        self.arg = arg


class Bdp(object):
    __slots__ = ['path','sid','mole_id','mole_struct','create_path']

    def __init__(self,path,sid):
        self.path = path
        self.create_path = os.path.dirname(os.path.dirname(path))
        self.sid = sid
        self.mole_id = []
        self.mole_struct = []
        # 初始化数据
        self.read_formula()
        self.init_path()
        
    @property
    def chain_num(self):
        s = self.read_struct()
        first_model = s[0]
        return len(first_model)

    @property
    def header_infos(self):
        return parse_pdb_header(self.path)

    @property
    def doi(self):
        _journal = self.header_infos['journal']
        _doi = _journal.split()[-1]
        if _doi.startswith('10.'):
            return _doi
        else:
            raise DoiError(f'current pdb does not doi, {_journal}')

    def mkdir(self,pname,md = False):
        """mkdir [summary]
        
        [extended_summary]
        
        :param pname: [description]
        :type pname: [type]
        :param md: [dict], defaults to False
        :type md: [dict], optional
        """
        def wrapper(func):
            def deco(self,*arg,**kwargs):
                create_path = os.path.join(self.create_path,pname)
                if os.path.exists(create_path):
                    ...
                else:
                    os.mkdir(create_path)
                # create markdown
                mdc = mdcheck(md)
                header_info = f'''---
    title: {mdc['title']}
    date: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}
    type: {mdc['typeinfo']}
    ---
    [TOC]
    # Description
    {mdc['description']}
    '''
                create_path_md = os.path.join(create_path,'README.md')
                with open(create_path_md,'w+') as f:
                    f.write(header_info)
                func(self,*arg,**kwargs)
            return deco
        return wrapper

    def value_check(self):
        try:
            if len(self.mole_id) == len(self.mole_struct):
                print('check pass!')
            else:
                raise ValueError('molecule identifier do not match formula identifier! manual check')
        except :
            raise ValueError

    @staticmethod
    def read_line_list(first_column,path):
        """read_line_list 读取第一列为first_column字符串的行，存为列表返回

        [extended_summary]

        Arguments:
            first_column {[string]} -- [pdb文件内容第一列常常为大写]
            path {[string]} -- [路径]

        Returns:
            [iter] -- [含有first_column字符串的生成器]
        """
        stringiter = Bdp.stringlinesiter(file = path)
        stringiter = map(lambda x:x.strip(), stringiter)
        return filter(lambda x:x.split()[0]==first_column,stringiter)

    @staticmethod
    def stringlinesiter(file):
        with open(file, 'r+', encoding='utf-8') as f:
            yield from f.readlines()

    @mkdir(self = 'self',pname = 'split_molecule')
    def split_molecule(self,residue):
        """split_molecule [split molecule or residue line save to file]
        
        [extended_summary]
        
        :param residue: [select residue]
        :type residue: [string]
        :extract_chain: extract molecule corresponding chain
        :type extract_chain: bool
        """
        remove_list,residue = [],str(residue)
        s = self.read_struct()
        chainlist = self.getchain_str()
        # 对每条链上都与这个小分子结合的链上的结合信息都提取出来
        for i in chainlist:
            io = PDBIO()
            io.set_structure(s[0][i])
            sidc = self.sid.replace('.pdb', '') if '.pdb' in self.sid else self.sid
            savename = f'{sidc}_{i}_{residue}.pdb'
            path = os.path.join(self.create_path,'split_molecule',savename)
            io.save(path,ResidueSelect(residue))
            if os.path.exists(path):
                with open(path,'r+') as f:
                    content = f.readline()
                    if content.strip() == 'END':
                        # 该链没有小分子，删除文件
                        # print(f'{savename} this chain have not molecule remove it')
                        remove_list.append(path)
        for i in remove_list:
            os.remove(i) # remove the empty molecule file


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

    def site_exists(self,residue,num):
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
                print(e,f'chain {i} not exist {residue}-{num}')
        return _l

    def _split_chain_save_struct(self,obj,sid):
        io = PDBIO()
        io.set_structure(obj)
        sidc = sid.replace('.pdb', '') if '.pdb' in sid else sid
        name,path = f'{sidc}_{obj.id}.pdb', self.create_path
        io.save(os.path.join(path,'split_chain',name),ChainSelect(obj.id))

    @mkdir(self = 'self',pname ='split_chain')
    def split_chain(self,extract='auto'):
        """split_chain [summary]
        
        [extended_summary]
        
        :param extract: [extract the chain], defaults to 'auto' extract all chains
        :type extract: [string], optional
        :raises ValueError: [save to pdb in current path in ./splitchain]
        """
        # 将传入的蛋白文件按照蛋白链分割并按照 pdbid_A_m.pdb 重新命名 A代表链，m表示该链含有小分子
        p = PDBParser(PERMISSIVE=1)
        s = p.get_structure(self.sid, self.path)
        if extract == 'auto':
            chainlist = self.getchain_str()
            for j in chainlist:
                self._split_chain_save_struct(obj=s[0][j],sid=self.sid)
        else:
            chainlist = self.getchain_str()
            if extract not in chainlist:
                raise ChainError(f'The {extract} chain not exists !')
            if extract in chainlist:
                extract = extract.upper()
                self._split_chain_save_struct(obj=s[0][extract],sid=self.sid)
            else:
                raise ValueError

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

    def get_chain(self,chain_str):
        s = self.read_struct()
        return s[0][chain_str]

    def init_path(self):
        _p = return_path(self.path)
        self.path = _p
        if os.path.isfile(_p) == False:
            raise FileNotFoundError(f'file {_p} not found')

    def read_formula(self):
        """read_formula read pdb file molecule information and identifier in this pdb file(just like residue)
        
        [extended_summary]
        
        :raises OSError: Do not find this file
        """
        _path = self.path
        try:
            with open(_path,'r') as f:
                text_line = f.readlines()
                _l = []
                for i in text_line:
                    li = i.split()
                    if li[0] == 'FORMUL':
                        _l.append(li)
                        self.mole_id.append(li[2])
                for ii in _l:
                    self.mole_struct.append(''.join(ii[3:]))
        except Exception as e:
            raise e

        def clean_pdb(self,outfile=None,path=None):
            __lst = Bdp.read_line_list('ATOM',self.path)
            _p = Path(path) if path else Path(self.path).absolute().parent
            _f = outfile if outfile else f"{self.sid}.clean.pdb"
            write_path = _p.joinpath(_f)
            if write_path.exists(): write_path.unlink() 
            with open(write_path,'w+') as f:
                for i in __lst:
                    f.write(f"{i}\n")

def renamepdbfile(path):
    """renamepdbfile 自动对path路径中批量下载的pdb中文件的文件名进行重命名

    [extended_summary]

    Arguments:
        path {[type]} -- [description]
    """
    _flist = get_path_contents(path)[0]
    for i in _flist:
        i = i.strip()
        basetuple = os.path.split(i)
        basedir = basetuple[0]
        _fn = basetuple[1]
        res = re.match(pattern='pdb([A-Za-z0-9]*).ent',string=_fn)
        # nfn = res.group(1)
        # newname =os.path.join(basedir,str(nfn)+'.pbd')
        # os.rename(src = i,dst = newname)
        try:
            nfn = res.group(1)
            newname =os.path.join(basedir,str(nfn)+'.pdb')
            os.rename(src = i,dst = newname)
        except IndexError:
            raise IndexError(f'Could not find pdbid in this filename: {i}')
        except AttributeError:
            continue
        except FileNotFoundError as e:
            if os.path.exists(i):
                raise e
            elif os.path.exists(newname):
                continue
            else:
                raise e



def pdbupdate(path):
    pl = PDBList(pdb=path)
    pl.update_pdb()


