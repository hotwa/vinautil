#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :pymol_api
@Description:       : some pymol api remake
@Date               :2023/1/1 15:06:35
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
import sys
from pathlib import Path
from typing import Optional, Tuple, NoReturn, Dict, List, Iterable, Callable, Any
from pymol import cmd, finish_launching, rpc
from utils.typecheck import typeassert

"""
cmd.centerofmass('organic and chain A') # 重心
# finish_launching(['pymol','-R']) # without GUI pymol.finish_launching(['pymol', '-cq']) # for windows
"""

class api():
    std_types = ['CYS', 'ILE', 'SER', 'VAL', 'GLN', 'LYS', 'ASN', 'PRO', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP',
                 'LEU', 'ARG', 'TRP', 'GLU', 'TYR', 'MET', 'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
    regular_residue = [i.upper() for i in 'Gly,Ala,Val,Leu,Ile,Pro,Phe,Tyr,Trp,Ser,Thr,Cys,Met,Asn,Gln,Asp,Glu,Lys,Arg,His'.split(',')]
    """注意: 在PyMOL中，"organic"和"hetatm"是有区别的，"organic"只包含碳氢化合物，而"hetatm"包含所有非蛋白质 / 核酸 / 多肽分子。
    "organic": 找到含碳的小分子，但不匹配已知的聚合物。使用方法：organic。例如：indicate organic。
    "all": 当前加载到 PyMOL 中的所有原子。
    "none": 空选择。
    "enabled": 启用对象中的原子。
    "sele": 命名选择或对象 "sele"，但仅当它与其他操作符的名称不冲突时。
    "%sele": 命名选择或对象 "sele"。推荐使用，以避免模糊性。
    atom: selects all atoms
    chain: selects atoms based on their chain identifier
    element: selects atoms based on their element type (e.g. "C" for carbon, "H" for hydrogen)
    resn: This selection name selects atoms based on their residue name (e.g. "resn ALA" selects all atoms in alanine residues).
    resi: selects atoms based on their residue number
    name: selects atoms based on their name (e.g. "CA" for alpha carbon)
    hetatm: selects atoms that are part of HETATM records in a PDB file
    polymer: selects atoms that are part of known polymers
    solvent: selects atoms that are part of solvent molecules
    organic: selects atoms that are part of carbon-containing molecules that do not match known polymers
    inorganic: selects atoms that are not part of carbon-containing molecules or known polymers
    br.(byres) all: This selection name selects all atoms that are within a certain distance of another atom or group of atoms (e.g. "br. all within 5 of organic" selects all atoms within 5 angstroms of any organic small molecules).
    """

    @staticmethod
    def pymol_start(args: list, rpc:bool = False, **kwargs) -> NoReturn:
        """start pymol
        pymol -q: 启动无GUI的pymol，同时忽略pymolrc文件（如果存在）。
        pymol -cq: 启动无GUI的pymol，并忽略pymolrc文件（如果存在）。这是在windows上启动无GUI的pymol的常用方法。
        finish_launching(['pymol','-cq']) # for windows start pymol without GUI
        finish_launching(['pymol','-c']) # for linux start pymol with GUI
        MacOS: not support finish_launching
        use subprocess.Popen(['pymol','-cq', '-R']) instead open XMLRPC server for interaction
        default port: 9123 # not support change port
        """
        start_cmd = {
            'win32': finish_launching,
            'linux': finish_launching,
            'darwin': 'MacOS not support finish_launching, use subprocess.Popen instead',
        }
        if isinstance(start_cmd[sys.platform], Callable):
            start_cmd[sys.platform](args)
        elif isinstance(start_cmd[sys.platform], str):
            print(start_cmd[sys.platform])
        else:
            raise Exception(f'Not support {sys.platform} platform!')
        hostname = kwargs.get('hostname') if kwargs.get('hostname') else 'localhost'
        _xmlPort = kwargs.get('port') if kwargs.get('port') else 9123
        _nPortsToTry = kwargs.get('nPortsToTry') if kwargs.get('nPortsToTry') else 5
        if rpc:
            rpc.launch_XMLRPC(hostname, _xmlPort, _nPortsToTry)

    @staticmethod
    def pymol_exit() -> NoReturn:
        '''
        @Description:       : exit pymol
        :return:           : NoReturn
        '''
        cmd.quit()

    def load_file(self, file: Path = None, clean:bool = True) -> NoReturn:
        cmd.delete("all")
        file = file if file else self.file
        cmd.load(file.as_posix())
        if clean: cmd.remove('solvent metals')  # remove solvent and metals
    def save_molecule(self, ligandid: str, chain:str, format: str = "mol2", file:Path = None) -> NoReturn:
        ligandid, chain = ligandid.upper(), chain.upper()
        cmd.delete("all")
        cmd.select(f"{ligandid}", f"chain {chain} and resn {ligandid}")
        if file is None:
            file = self.file.parent.joinpath(f'{ligandid}_{chain}.{format}').as_posix()
        else:
            file = file.as_posix()
        cmd.save(file, ligandid)
        cmd.delete("all") # clean pymol window

    @staticmethod
    def AddPolarHydrogens(selection:str, out_file: Path = None, fmt:str = 'mol2') -> Optional[str]:
        """
        The file format is automatically chosen if the extesion is one of
        the supported output formats: pdb, pqr, mol, sdf, pkl, pkla, mmd, out,
        dat, mmod, cif, pov, png, pse, psw, aln, fasta, obj, mtl, wrl, dae, idtf,
        or mol2.
        :param selection: pymol grammar selector
        :param out_file: output file
        :return:
        openbabel method:
        omol = list(pybel.readfile(format = 'mol2', filename = file))[0]
        omol.OBMol.AddPolarHydrogens()
        omol.write('mol2',out_file,overwrite=True)
        """
        cmd.delete("all")
        cmd.h_add(selection) # add all hydrogens in this molecular
        cmd.remove(f"{selection} & hydro & not nbr. (don.|acc.)") # remove no-polar-hydrogen
        if out_file:
            cmd.save(filename=out_file.as_posix(), selection=selection)
        else:
            return cmd.get_str(format=fmt, selection=selection)
    @staticmethod
    def pymol_get_sequence(resn:str, chain:str) -> List[str]:
        """
            获取链{chain}残基名称为{resn}的残基序列。
        """
        res = []
        cmd.iterate(selection=f'resn {resn} and chain {chain}', expression=lambda atom: res.append(atom.resi))
        res_seq = list(set(res))
        if len(res_seq) > 1: # 在同一条链上可以有两个以上相同的小分子
            print(f"Warning: resn {resn} and chain {chain} has more than one residue sequence: output {res_seq}")
        return res_seq

    @staticmethod
    def pymol_get_around(resn:str, chain:str, distance:float, extend: bool = True, add_self:bool = True) -> str:
        """
            获取链{chain}残基名称为{resn}的周围{distance}A的pdb字符。
            :param resn: 残基名称
            :param chain: 链
            :param distance: 距离
            :param add_self: 是否添加resn选择对象自身
            :return: pdb字符
file = 'pdb7axn.ent.gz'

resn = 'S6B'
chain = 'A'
distance = 3.00
cmd.delete("all")
cmd.load(file.as_posix())
cmd.remove('metal solvent')
cmd.select('sele', f'resn {resn} and chain {chain}')
cmd.center('sele')
cmd.orient('sele')
obj_name = 'around_and_self'
# cmd.create(obj_name, f'(resn {resn} and chain {chain}) around {distance}') # 不扩展残基
# cmd.create(obj_name, f'byres (resn {resn} and chain {chain}) around {distance}') # 扩展残基
# cmd.create(obj_name, f'(resn {resn} and chain {chain}) around {distance} or (resn {resn} and chain {chain})') # 选择resn的对象和resn对象周围3A的原子
cmd.create(obj_name, f'byres (resn {resn} and chain {chain}) around {distance} or (resn {resn} and chain {chain})') # 选择resn的对象和resn对象周围3A的原子扩展至残基的原子
cmd.hide('everything', 'all')
cmd.show('lines', obj_name)
pdbstr = cmd.get_pdbstr(selection=obj_name)
print(pdbstr)
        """
        obj_name = f'chain_{chain}_{resn}_around_{round(distance,2)}'
        if add_self:
            if extend:
                cmd.create(obj_name,
                           f'byres (resn {resn} and chain {chain}) around {distance} or (resn {resn} and chain {chain})')
            else:
                cmd.create(obj_name, f'(resn {resn} and chain {chain}) around {distance} or (resn {resn} and chain {chain})')
        else:
            if extend:
                cmd.create(obj_name, f'byres (resn {resn} and chain {chain}) around {distance}')
            else:
                cmd.create(obj_name, f'(resn {resn} and chain {chain}) around {distance}')
        pdb_str = cmd.get_pdbstr(selection=obj_name)
        cmd.delete(obj_name) # delete the object
        return pdb_str
    @staticmethod
    def pymol_get_centerofmass(object_name) -> List[float]:
        """
            获取 pymol 中对象的重心。

            这是一个用来获取 pymol 中对象的重心的函数坐标

            :param object_name: 对象的名称
            :type object_name: str，optional
            :return: 重心的坐标。
        """
        return cmd.centerofmass(object_name)

    @staticmethod
    def pymol_print_exec(object_name: str = 'organic', atom_attribute: str = 'resn') -> List[str]:
        """
            打印 pymol 中对象的属性。

            这是一个用来打印 pymol 中对象的属性的函数。

            :param object_name: 对象的名称，默认为 'organic'。
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

    @staticmethod
    def get_coords(selection: str):
        return cmd.get_coords(selection)
    @staticmethod
    def pymol_get_resn(object_name: str = 'all') -> Tuple[str]:
        resn_list = []
        cmd.iterate(selection=object_name, expression=lambda atom: resn_list.append(atom.resn))
        return tuple(set(resn_list))

    @staticmethod
    def Mutagenesis_site(filename: Path, mutation_type: str, site: int, outfile: Path = None) -> Path:
        """Mutagenesis_site Mutagenesis residue in site
            residue site have mutil conformations, need to select one conformations, some error accured.
            pass
        _extended_summary_

        Arguments:
            filename {str} -- PDB file format
            mutation_type {str} -- 'VAL' for ALA TO VAL; 'ALA' for any/not ALA to ALA; 'GLY' for ALA to GLY
            site {int} -- residue site in pdbfile

        Keyword Arguments:
            outfile {str} -- _description_ (default: {None})

        Raises:
            ValueError: not one object in PDBs,need to fix

        Returns:
            str -- save mutagenesis file path
        """
        p = Path(filename)
        savename = p.stem + f"_{site}_mutation.pdb"
        _out_file = Path(outfile) if outfile else p.absolute().parent.joinpath(savename)
        if not _out_file.absolute().parent.exists(): _out_file.absolute().parent.mkdir(parents=True)
        cmd.delete('all')  # ! clean up
        cmd.load(filename)
        PDBs = cmd.get_names()
        # Get the ID numbers of c-alpha (CA) atoms of all residues
        if len(PDBs) == 1:
            PDB = PDBs[0]
        else:
            raise ValueError(f'this pdb have more than one object!PDBs:{PDBs}')
        CAindex = cmd.identify(f"{PDB} and name CA")
        pdbstrList = [cmd.get_pdbstr("%s and id %s" % (PDB, CAid)).splitlines() for CAid in CAindex]
        ProtChainResiList = [[PDB, i[0][21], i[0][22:26].strip()] for i in pdbstrList]
        for i, j, k in ProtChainResiList:
            if int(k) == int(site):
                cmd.wizard("mutagenesis")
                # print(i,j,k)
                cmd.refresh_wizard()
                cmd.get_wizard().set_mode(mutation_type)
                ##Possible mutation_type could be:
                ##'VAL' for ALA TO VAL
                ##'ALA' for any/not ALA to ALA
                ##'GLY' for ALA to GLY
                # 'selection' will select each residue that you have selected
                # on ProtChainResiList above using the PDBid,chain,and residue
                # present on your pdb file.If you didn't select a range on
                # ProteinChainResiList, it will do the mutation on all the residues
                # present in your protein.
                selection = f"/{i}//{j}/{k}"
                # Print selection to check the output
                # print(selection)
                # Selects where to place the mutation
                cmd.get_wizard().do_select(selection)
                ##Applies mutation
                cmd.get_wizard().apply()
        # Save each mutation and reinitialize the session before the next mutation
        ##to have pdb files only with the residue-specific single-point mutation you were interested.
        cmd.set_wizard("done")
        cmd.save(_out_file.as_posix(), f"{PDB}")
        cmd.delete('all')  # Reinitialize PyMOL to default settings.
        return _out_file