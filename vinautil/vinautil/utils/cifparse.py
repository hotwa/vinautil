from pathlib import Path
import gzip
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Union, Optional, NoReturn, Generator, Iterable, Callable, Any
from Bio.PDB.MMCIF2Dict import MMCIF2Dict as mmcifdict
from pdbecif.mmcif_tools import MMCIF2Dict

try:
    from .typecheck import typeassert
    from .pdbparse_module import Bdp
except (ImportError, ModuleNotFoundError, AttributeError, NameError, SyntaxError):
    from typecheck import typeassert
    from pdbparse_module import Bdp

@typeassert(file = Path)
@dataclass
class cif_parse():
    """
    parse cif file
    """
    file: Path
    pid: str = field(default_factory=lambda: None)
    file_content: str = field(init=False)
    cif_dict: dict = field(init=False)
    bio_cif_dict: Dict = field(init=False)
    pdbecif_dict: Dict = field(init=False)

    def __post_init__(self):
        if self.file.suffix == '.gz':
            self.file_content = self.__get_cif_gz_content()
        self.pdbecif_dict = MMCIF2Dict().parse(self.file.as_posix())
        self.bio_cif_dict = mmcifdict(self.file.as_posix())
        self.cif_dict = self.pdbecif_dict
        if self.pid is None:
            self.pid = self.__get_pid()

    def __get_cif_gz_content(self) -> str:
        if '.gz' in self.file.suffixes:
            with gzip.open(self.file.as_posix(), 'rt') as f:
                content = f.read()
        else:
            raise ValueError(f"{self.file.as_posix()} is not a gzip file, or not suffix with '.gz'")
        if '.cif' in Path(self.file.stem).suffix:
            rename = self.file.stem
        else:
            rename = self.file.name.split('.')[0] + '.cif'
        ungzip_file = self.file.parent.joinpath(rename).as_posix()
        with open(ungzip_file, 'w', encoding='utf-8') as f:
            f.write(content)
            self.file = Path(ungzip_file)
        return content

    def search_keylike(self, key: str) -> List[str]:
        return [{i:self.cif_dict[i]} for i in self.cif_dict.keys() if key in i]

    def __get_pid(self) -> str:
        if len(self.bio_cif_dict['_entry.id']) == 1:
            return self.bio_cif_dict['_entry.id'][0]
        else:
            raise ValueError(f"More than one PDB ID in {self.file.as_posix()}")

    @property
    def conn_infos(self) -> Dict:
        return self.cif_dict[self.pid]['_struct_conn']

    @property
    def auth_comp_id(self) -> List[str]:
        return list(set(self.bio_cif_dict['_struct_site.pdbx_auth_comp_id']))

    @property
    def covalent_bonds(self) -> bool:
        if 'covale' in self.bio_cif_dict['_struct_conn.conn_type_id']:
            return True
        else:
            return False








"""
_pdbx_validate_symm_contact.id
_pdbx_validate_torsion.auth_asym_id
_pdbx_struct_assembly_gen.asym_id_list
_struct_conn.conn_type_id
_struct_site.pdbx_auth_asym_id
_struct_site.pdbx_auth_comp_id
_struct_site.details
_struct_conn.conn_type_id字段还有许多其他连接类型，包括：

covale: 共价键
disulf: 二硫键
covale1: 共价键1
covale2: 共价键2
covale3: 共价键3
hbond: 电子转移键
metalc: 金属配位键
metalc1: 金属配位键1
metalc2: 金属配位键2

disulf: 双硫键连接。
covale: 共价连接。
metalc: 金属连接。
saltbr: 离子之间的相互作用。
hydrogen: 氢键连接。
water: 水分子之间的连接。
helix: 螺旋结构。
turn: 旋转结构。
cis: 同侧螺旋结构。
trans: 异侧螺旋结构。
xlink: 蛋白质间的交联。
covalent: 共价连接。
"""