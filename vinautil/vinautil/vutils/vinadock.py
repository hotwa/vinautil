from pathlib import Path
from vina import Vina

# ------dockvina1.2.3 example script------

def dockvina(receptor='1iep_receptor.pdbqt',ligand='1iep_ligand.pdbqt',center=[15.190, 53.903, 16.917],box_size=[20, 20, 20],exhaustiveness=32,n_poses=20,out_n_poses = 5):
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
    v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)
    # Dock the ligand
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=out_n_poses, overwrite=True)