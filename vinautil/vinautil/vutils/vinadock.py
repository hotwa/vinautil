from pathlib import Path
from vina import Vina

def dockvina(receptor, ligand, center, box_size, exhaustiveness=32,n_poses=20,out_n_poses = 5):
    out_stem = {Path(receptor).stem}--{Path(ligand).stem}
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
    v.write_pose(f'{out_stem}_minimized.pdbqt', overwrite=True)
    # Dock the ligand
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(f'{out_stem}.pdbqt', n_poses=out_n_poses, overwrite=True)