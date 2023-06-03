import argparse
from pathlib import Path
import os,sys
import subprocess
import datetime
from loguru import logger
from vina import Vina
from pymol import cmd
from typing import List, float
from vinautil.pymolutils.mutagenesis import Mutagenesis_site

conda_prefix = Path(os.environ.get('CONDA_PREFIX'))
python2_interpreter = conda_prefix.joinpath('bin/python2')
python3_interpreter = conda_prefix.joinpath('bin/python3')
prepare_ligand4 = conda_prefix.joinpath('MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py')
prepare_receptor4 = conda_prefix.joinpath('MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py')
mk_prepare_ligand = conda_prefix.joinpath('bin/mk_prepare_ligand.py')

def dockvina(receptor:Path, ligand:Path, center:List[float], box_size:List[float], 
             exhaustiveness:int =32,n_poses:int =20,out_n_poses:int = 20):
    out_stem = {receptor.stem}--{ligand.stem}
    out_dir = receptor.parent
    v = Vina(sf_name='vina')
    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)
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
    v.write_poses(out_dir.joinpath(f'{out_stem}.pdbqt').as_posix(), n_poses=out_n_poses, overwrite=True)
    return v

def SCARdock(receptor, ligand, chain, site):
    # command line SCARdock
    print('''SCARdock method DOI: 10.1021/acs.jcim.6b00334
SCARdock screening server(https://scardock.com) DOI: 10.1021/acsomega.2c08147
lab site: http://liugroup.site
author: Lingyu Zeng mail: pylyzeng@gmail.com
''')
    # protein Mutagenesis (GLY)
    receptor = Path(receptor)
    if not receptor.exists():
        raise FileNotFoundError(f'{receptor} not found')
    muta_receptor = receptor.parent.joinpath(f'{receptor.stem}_{site}G.pdb')
    Mutagenesis_site(filename=receptor, mutation_type='GLY', site=int(site), outfile= muta_receptor)
    # prepare receptor dock file
    receptor_pdbqt = receptor.parent.joinpath(f"{muta_receptor.stem}.pdbqt")
    CMD_ = f'{python2_interpreter} {prepare_receptor4} -r {receptor.as_posix()} -o {receptor_pdbqt.as_posix()} -A checkhydrogens'
    p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
    while p.poll() is None:  # 表示进程还在运行
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
    # prepare ligand dock file
    ligand = Path(ligand)
    ligand_pdbqt = ligand.parent.joinpath(f"{ligand.stem}.pdbqt")
    if not ligand.exists():
        raise FileNotFoundError(f'{ligand} not found')
    CMD_ = f'{python3_interpreter} {mk_prepare_ligand} -i {ligand.as_posix()} -o {ligand_pdbqt.as_posix()}'
    p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
    while p.poll() is None:  # 表示进程还在运行
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
    # box center, covalent beta carbon coordinate
    cmd.reinitialize('everything')
    cmd.load(receptor.as_posix())
    coordinates = []
    cmd.select('res_covalent_atoms', f'chain {chain} and resi {site} and name CB') # select residue covalent site beta carbon, be careful GLY have no side chain with beta carbon
    cmd.iterate_state('0', 'res_covalent_atoms', 'coordinates.append([x,y,z])', space=locals())
    if len(coordinates) == 1:
        center = coordinates[0]
    else:
        raise ValueError(f'covalent beta carbon site {site} not found!')
    # SCARdock docking
    dockvina(receptor=receptor_pdbqt, ligand=ligand_pdbqt, center=center, box_size=[40, 40, 40], exhaustiveness=32,n_poses=20,out_n_poses = 20)

    
    
    
    
    
    
    
    
    
    
    
    
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

    
    
    
