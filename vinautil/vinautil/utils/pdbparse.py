from pathlib import Path
import gzip
import re
import time
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Union, Optional, NoReturn, Generator, Iterable, Callable, Any, Set
from enum import Enum
from myutils.pdbparse import Bdp
"""
prepare pdb file:
rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ /home/mmme/alist/local_storage/PDBdatabase/pdb/
rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ /home/mmme/alist/local_storage/PDBdatabase/cif/
transfer pdb file to pdbqt file:
ADFRsuite-1.0/bin/relpath_prepare_receptor4.bat -r test.scar.pdb -o test.scar.pdbqt -A checkhydrogens
"""

# here = Path(__file__).parent.resolve()
# import sys
# sys.path.insert(0, here)
from utils.typecheck import typeassert
from utils.pymol import api

class CovalentSourceType(Enum):
    PYMOL = 0 # use pymol show lines see covalent bond
    PAPER = 1 # from sci paper get covalent information
    LINK = 2 # from pdb file LINK record get covalent bond
    def __str__(self):
        return 'how to judge covalent bond is class CovalentType'
    def is_pymol(self) -> bool:
        return self == CovalentSourceType.PYMOL

    def is_paper(self) -> bool:
        return self == CovalentSourceType.PAPER

    def is_link(self) -> bool:
        return self == CovalentSourceType.LINK

class ConvalentType:
    SMALL_MOLECULE_SIDE = 'molecule_side'
    PROTEIN_SIDE = 'protein_side'
    def __init__(self, type: str):
        if type not in (self.SMALL_MOLECULE_SIDE, self.PROTEIN_SIDE):
            raise ValueError(f"Invalid ConvalentType: {type}")
        self.type = type

    def __repr__(self):
        return f"ConvalentType({self.type!r})"

    def __doc__(self):
        return f"ConvalentType representing a covalent bond on either the molecule side or the protein side."

    def __str__(self):
        return self.type

    @classmethod
    def is_small_molecule_side(cls, val: Union[str, 'ConvalentType']):
        return val == ConvalentType.SMALL_MOLECULE_SIDE or val == ConvalentType.SMALL_MOLECULE_SIDE

    @classmethod
    def is_protein_side(cls, val: Union[str, 'ConvalentType']):
        return val == ConvalentType.PROTEIN_SIDE or val == ConvalentType.PROTEIN_SIDE

@dataclass
class CovalentMsg():
    sourcetype: CovalentSourceType
    content: Optional[List[Dict]]
    """type: Union[ConvalentType, type(None)] = field(default=None)
"""
    def is_pymol(self) -> bool:
        return self.sourcetype == CovalentSourceType.PYMOL

    def is_paper(self) -> bool:
        return self.sourcetype == CovalentSourceType.PAPER

    def is_link(self) -> bool:
        return self.sourcetype == CovalentSourceType.LINK

@dataclass
class molecule():
    """
    A class representing a molecule.

    Attributes:
    resn (str): The residue name of the molecule.
    auth_label_chain (str): The author-assigned label and chain identifier of the molecule.
    seq (int): The sequence number of the molecule.
    cov (Union[CovalentMsg, None]): The covalent bond information of the molecule (if any). Default value is None.
    """
    resn: str
    auth_label_chain: str
    seq: int
    cov: Union[CovalentMsg, None] = field(default=None)


@typeassert(file=Path)
@dataclass
class pdb_parse(api):
    """
    Class for parsing PDB files.
    file: Path
        The file path of the PDB file to be parsed.
    pid: str
        The PDB ID of the molecule. This field is initialized to be empty and will be filled when the PDB file is parsed.
    file_content: List[str]
        A list of strings representing the lines of the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    hetatm: Tuple[str]
        A tuple of strings representing HETATM records in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    organic: Tuple[str]
        A tuple of strings representing organic HETATM records in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    metals: Tuple[str]
        A tuple of strings representing metal HETATM records in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    inorganic: Tuple[str]
        A tuple of strings representing inorganic HETATM records in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    chain: Tuple[str]
        A tuple of strings representing chain records in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    mole: List[molecule]
        A list of `molecule` objects representing the molecules in the PDB file. This field is initialized to be empty and will be filled when the PDB file is parsed.
    """
    file: Path
    pid: str = field(init=False)
    file_content: List[str] = field(init=False)
    hetatm: Tuple[str] = field(init=False)
    organic: Tuple[str] = field(init=False)
    metals: Tuple[str] = field(init=False)
    inorganic: Tuple[str] = field(init=False)
    chain: Tuple[str] = field(init=False)
    mole: List[molecule] = field(init=False)
    pymol_start_args: List[str] = field(default_factory=lambda: ["pymol", "-M", "-A1"])

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.file.as_posix()}')"
    def __post_init__(self):
        self.pid = self.match_pid(string=self.file.name)
        self.file_content = list(self.get_file_content())
        self.pymol_start(args=self.pymol_start_args)
        self.load_file() # load file in pymol
        self.__init_infos() # init base infos

    def __init_infos(self) -> NoReturn:
        self.chain = list(set(self.pymol_print_exec(atom_attribute='chain')))
        self.organic = self.pymol_get_resn(object_name='organic')
        self.hetatm = self.pymol_get_resn(object_name='hetatm')
        self.metals = self.pymol_get_resn(object_name='metals')
        self.inorganic = self.pymol_get_resn(object_name='inorganic')
        self.mole = self._get_molecule()

    def _get_molecule(self) -> List[molecule]:
        lst = []
        for c in self.chain:
            for m in self.organic:
                # print(f'in {self.file.name}, try get {c} chain {m} molecule')
                seq = self.pymol_get_sequence(resn=m, chain=c)
                if len(seq) == 1 : # one molecule is in this chain
                    lst.append(molecule(resn=m, auth_label_chain=c, seq=seq, cov=self._search_covalent_bond()))
                elif len(seq) > 1: # same chemical molecule have more than one molecule in this chain
                    for s in seq:
                        lst.append(molecule(resn=m, auth_label_chain=c, seq=s, cov=self._search_covalent_bond()))
                else:
                    pass # the molecule not in this chain
        return lst
    def _search_covalent_bond(self, distance:float = 3.00) -> List[Dict[str, str]]:
        """
        search covalent bond from pymol, around all molecule 3.00 ANG
        """
        # CovalentSourceType.PYMOL
        for c in self.chain:
            for mol in self.organic:
                pdb_str = self.pymol_get_around(resn=mol, chain=c, distance=distance)
                # print(pdb_str)
                search_cov = self.pdbstr_covalent(pdb_str=pdb_str)
                # print(search_cov)
        return search_cov

    @staticmethod
    def pdbstr_covalent(pdb_str: str) -> List[Dict[str, str]]:
        """
        get covalent bond from pdb string, bond connect ATOM line and HETATM line
        :param pdb_str:
        :return: list of dict
        """
        pdbstr_lines = pdb_str.splitlines()
        cov_bonds = []
        for line in pdbstr_lines:
            if line.startswith('CONECT'):
                conect_info = line
                if len(line.split()) == 3:
                    atoms = list(map(lambda x: int(x), line.split()[1:]))
                    atom_lines = list(
                        filter(lambda x: x.startswith('ATOM') or x.startswith('HETATM'), pdbstr_lines))
                    atoms_info = list(map(lambda x: x.strip(), atom_lines))
                    lst = [conect_info, atoms_info[atoms[0] - 1], atoms_info[atoms[1] - 1]]
                    cov_bond = {i.split()[0]:i for i in lst}
                    cov_bonds.append(cov_bond)
        return_value = [i for i in cov_bonds if len(i.keys()) == 3]
        return return_value

    def _get_file_content_iterable(self) -> Iterable[str]:
        with open(self.file.as_posix(), 'r') as f:
            for line in f:
                yield line

    def get_file_content(self) -> Iterable[str]:
        if '.ent' and '.gz' in self.file.suffixes:
            self._ungzip_pdb_gz()
            return self._get_file_content_iterable()
        elif '.pdb' in self.file.suffix:
            return self._get_file_content_iterable()
        elif '.ent' in self.file.suffix:
            return self._get_file_content_iterable()
        else:
            raise ValueError(f"File type error: {self.file}, only support .ent, .pdb or .ent.gz")

    def _ungzip_pdb_gz(self) -> NoReturn:
        fmt = 'pdb'
        with gzip.open(self.file.as_posix(), 'rt') as f_in:
            if self.pid is None:
                ungzip_file_name = f"unkown_PDB_ID_{int(time.time())}.{fmt}" # 如果达到秒级运算，则需要修改删除int，否则会覆盖文件
            else:
                ungzip_file_name = f"{self.pid}.{fmt}"
            ungzip_file = self.file.parent.joinpath(ungzip_file_name)
            with open(ungzip_file.as_posix(), 'w', encoding='utf-8') as f_out:
                for line in f_in:
                    f_out.write(line)
        self.file = ungzip_file

    def _get_link_info(self) -> List[str]:
        return self._get_start_with("LINK")

    def _get_start_with(self, start_word: str) -> List[str]:
        return [line for line in self.file_content if line.startswith(start_word)]

    @staticmethod
    def match_pid(string) -> Optional[str]:
        # pattern = r'(pdb)?(?P<pdbid>[A-Za-z\d]{4}).(ent|pdb)(.gz)?'
        # ['r3w2.pdb', 'r3w2.pdb.gz', 'pdbr3w2.ent', 'pdbr3w2.ent.gz', 'pdbr3w2'] # pdbid = r3w2 test
        pattern = r'(pdb)?(?P<pdbid>[A-Za-z\d]{4})(\.(ent|pdb)(\.gz)?)?'
        match = re.match(pattern, string)
        if match:
            return match.group('pdbid')
        return None

    def link_info(self):
        pass


def get_gz_content(file: Path, ungzip:bool = True) -> str:
    if '.gz' in file.suffixes:
        with gzip.open(file.as_posix(), 'rt') as f:
            content = f.read()
    else:
        raise ValueError(f"{file.as_posix()} is not a gzip file, or not suffix with '.gz'")
    if ungzip:
        rename = file.name.split('.')[0]
        ungzip_file = file.parent.joinpath(rename).as_posix()
        with open(ungzip_file, 'w', encoding='utf-8') as f:
            f.write(content)
    return content




"""
共价键记录不一定只在LINK开头的信息中，还有可能在非LINK开头的信息中，例如：REMARK 500
PDB ID: 3sh8
REMARK 500 THE FOLLOWING ATOMS ARE IN CLOSE CONTACT.                            
REMARK 500                                                                      
REMARK 500  ATM1  RES C  SSEQI   ATM2  RES C  SSEQI           DISTANCE          
REMARK 500   OG   SER A    70     C8   CED A     1              1.90            
REMARK 500   OG   SER B    70     C8   CED B     1              2.03            
REMARK 500                                                                      
REMARK 500 REMARK: NULL                                                         
REMARK 500                                                                      
REMARK 500 GEOMETRY AND STEREOCHEMISTRY                                         
REMARK 500 SUBTOPIC: COVALENT BOND LENGTHS                                      
"""