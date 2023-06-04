import argparse
from pathlib import Path
import os,sys
import subprocess
import datetime
from loguru import logger

from pymol import cmd
from openbabel import pybel
from typing import List
from vinautil.vutils.obabel import PDBQTtoMol2, PDBQTparser
from vinautil.vutils.spyrmsd_load import symmrmsd_mol2_list
from vinautil.vina import Vina
from vinautil.pymolutils.mutagenesis import Mutagenesis_site

conda_prefix = Path(os.environ.get('CONDA_PREFIX'))
python2_interpreter = conda_prefix.joinpath('bin/python2')
python3_interpreter = conda_prefix.joinpath('bin/python3')
prepare_ligand4 = conda_prefix.joinpath('MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py')
prepare_receptor4 = conda_prefix.joinpath('MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py')
mk_prepare_ligand = conda_prefix.joinpath('bin/mk_prepare_ligand.py')

def dockvina(receptor:Path, ligand:Path, center:List[float], box_size:List[float], 
             exhaustiveness:int =32,n_poses:int =20,out_n_poses:int = 20):
    out_stem = f'{receptor.stem}--{ligand.stem}'
    out_dir = receptor.parent
    v = Vina(sf_name='vina')
    v.set_receptor(receptor.as_posix())
    v.set_ligand_from_file(ligand.as_posix())
    v.compute_vina_maps(center=center, box_size=box_size)
    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])
    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(out_dir.joinpath(f'{out_stem}_minimized.pdbqt').as_posix(), overwrite=True)
    # Dock the ligand
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    docked_file = out_dir.joinpath(f'{out_stem}.pdbqt')
    v.write_poses(docked_file.as_posix(), n_poses=out_n_poses, overwrite=True)
    return docked_file

def SCARdock():
    # cmd lineparser
    parser = argparse.ArgumentParser(description='SCARdock instance')
    parser.add_argument('-r', '--receptor',  metavar="recepotr file", nargs='?', default=sys.stdin, help='recepotr file, support pdb')
    parser.add_argument('-l', '--ligand',  metavar="ligand file", nargs='?', default=sys.stdin, help='ligand file, molecule file (MOL2, SDF,...)(use meeko prepare)')
    parser.add_argument('-s', '--site',  metavar="residue covalent site", nargs='?', default=sys.stdin, help='residue covalent site')
    parser.add_argument('-c', '--chain',  metavar="covalent chain ID", nargs='?', default=sys.stdin, help='covalent chain ID')
    parser.add_argument('-log', '--log_dir', metavar="output log directory", nargs='?', default='./', help='Relative Path')
    args = parser.parse_args()
    # 使用logger.add()方法，指定日志文件的路径和名称，以及编码方式
    log_file = Path(args.log_dir) / f'scardock.log'
    logger.add(log_file.as_posix(), encoding='utf-8')
    logger.info(f'SCARdock docking...')
    # command line SCARdock
    print('''SCARdock method DOI: 10.1021/acs.jcim.6b00334
SCARdock screening server(https://scardock.com) DOI: 10.1021/acsomega.2c08147
lab site: http://liugroup.site
author: Lingyu Zeng mail: pylyzeng@gmail.com
''')
    # clean pdb file
    receptor = cleanATOM(args.receptor) # same pyrosetta.toolbox cleanATOM
    if not receptor.exists():
        raise FileNotFoundError(f'{receptor.as_posix()} not found')
    # protein Mutagenesis (GLY)
    muta_receptor = receptor.parent.joinpath(f'{receptor.stem}_{args.site}G.pdb')
    Mutagenesis_site(filename=receptor, mutation_type='GLY', site=int(args.site), outfile= muta_receptor)
    # prepare receptor dock file
    receptor_pdbqt = receptor.parent.joinpath(f"{muta_receptor.stem}.pdbqt")
    CMD_ = f'{python2_interpreter} {prepare_receptor4} -r {receptor.as_posix()} -o {receptor_pdbqt.as_posix()} -A checkhydrogens'
    p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
    while p.poll() is None:  # 表示进程还在运行
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
    # prepare ligand dock file(mol2 file)
    ligand = Path(args.ligand)
    ligand_pdbqt = ligand.parent.joinpath(f"{ligand.stem}.pdbqt")
    if not ligand.exists():
        raise FileNotFoundError(f'{ligand} not found')
    # use openbabel add polar hydrogens
    molH = pybel.readfile('mol2', ligand.as_posix())
    molH = next(molH)
    molH.OBMol.DeleteHydrogens()
    molH.OBMol.AddPolarHydrogens()
    molH.write('mol2', ligand.as_posix(), overwrite=True)
    # use meeko backend prepare ligand, optional MGLtools prepare_ligand4.py
    CMD_ = f'{python3_interpreter} {mk_prepare_ligand} -i {ligand.as_posix()} -o {ligand_pdbqt.as_posix()}'
    p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
    while p.poll() is None:  # 表示进程还在运行
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
    # box center, covalent beta carbon coordinate
    cmd.reinitialize('everything')
    cmd.load(receptor.as_posix())
    coordinates = []
    cmd.select('res_covalent_atoms', f'chain {args.chain} and resi {args.site} and name CB') # select residue covalent site beta carbon, be careful GLY have no side chain with beta carbon
    cmd.iterate_state('0', 'res_covalent_atoms', 'coordinates.append([x,y,z])', space=locals())
    if len(coordinates) == 1:
        center = coordinates[0]
    else:
        raise ValueError(f'covalent beta carbon site {args.site} not found!')
    # SCARdock docking
    # ! be careful, orginal molecule coordinate should in docking box
    docked_file = dockvina(receptor=receptor_pdbqt, ligand=ligand_pdbqt, center=center, box_size=[40, 40, 40], exhaustiveness=32,n_poses=20,out_n_poses = 20)
    # restore mol2
    ins = PDBQTtoMol2(ligand.read_text(),ligand_pdbqt.read_text(),PDBQTparser(docked_file).get_modules())
    res_mol2 = ins.to_string() # mol2 string list
    # write mol2
    for n,m in enumerate(res_mol2):
        mol2_file = ligand.parent.joinpath(f'{ligand.stem}+{str(n+1).zfill(3)}.mol2')
        mol2_file.write_text(m)
    # cal RMSD(spyrmsd)
    res = symmrmsd_mol2_list(mol2_docked=res_mol2, mol2_ref=ligand.read_text())
    res_list = ['{:.2f}'.format(i) for i in res]
    print(f'''
SCARdock docking result: {docked_file}
RMSD: {res_list}''')


    
def cleanATOM(pdb_file, out_file=None, ext="_clean.pdb")->Path:
    """Extract all ATOM and TER records in a PDB file and write them to a new file.

    Args:
        pdb_file (str): Path of the PDB file from which ATOM and TER records
            will be extracted
        out_file (str): Optional argument to specify a particular output filename.
            Defaults to <pdb_file>.clean.pdb.
        ext (str): File extension to use for output file. Defaults to ".clean.pdb"
    """
    # find all ATOM and TER lines
    with open(pdb_file, "r") as fid:
        good = [l for l in fid if l.startswith(("ATOM", "TER"))]

    # default output file to <pdb_file>_clean.pdb
    if out_file is None:
        out_file = os.path.splitext(pdb_file)[0] + ext

    # write the selected records to a new file
    with open(out_file, "w") as fid:
        fid.writelines(good)
    return Path(out_file)
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SCARdock instance')
    parser.add_argument('-r', '--receptor',  metavar="recepotr file", nargs='?', default=sys.stdin, help='recepotr file, support pdb')
    parser.add_argument('-l', '--ligand',  metavar="ligand file", nargs='?', default=sys.stdin, help='ligand file, molecule file (MOL2, SDF,...)(use meeko prepare)')
    parser.add_argument('-s', '--site',  metavar="residue covalent site", nargs='?', default=sys.stdin, help='residue covalent site')
    parser.add_argument('-c', '--chain',  metavar="covalent chain ID", nargs='?', default=sys.stdin, help='covalent chain ID')
    parser.add_argument('-log', '--log_dir', metavar="output log directory", nargs='?', default='./', help='Relative Path')
    args = parser.parse_args()
    # 使用logger.add()方法，指定日志文件的路径和名称，以及编码方式
    log_file = Path(args.log_dir) / f'dock.log'
    logger.add(log_file.as_posix(), encoding='utf-8')
    logger.info(f'SCARdock docking...')
    SCARdock(receptor=args.receptor, ligand=args.ligand, chain=args.chain, site=args.site)

    
    
    
