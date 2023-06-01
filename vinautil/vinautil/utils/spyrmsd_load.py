
from spyrmsd.molecule import Molecule
from spyrmsd import io, rmsd
from spyrmsd.exceptions import NonIsomorphicGraphs
from numpy import nan as np_nan


def symmrmsd_mol2(mol2_ref: Molecule, mol2_docked: Molecule, minimize=False) -> float:
    ref = mol2_ref
    mol = mol2_docked
    ref.strip()  # remove Hydrogen atom
    mol.strip()  # remove Hydrogen atom
    coords_ref = ref.coordinates
    anum_ref = ref.atomicnums
    adj_ref = ref.adjacency_matrix
    coords = [mol.coordinates]
    anum = mol.atomicnums
    adj = mol.adjacency_matrix
    try:
        RMSD = rmsd.symmrmsd(
            coords_ref,
            coords,
            anum_ref,
            anum,
            adj_ref,
            adj,
            minimize=minimize
        )
        return RMSD[0]
    except NonIsomorphicGraphs:
        return np_nan

def rmsd_mol2(mol2_ref: Molecule, mol2_docked: Molecule, minimize=False) -> float:
    mol_ref = mol2_ref
    mol_docked = mol2_docked
    mol_ref.strip()  # remove Hydrogen atom
    mol_docked.strip()  # remove Hydrogen atom
    RMSD = rmsd.rmsd(
        coords1=mol_ref.coordinates,
        coords2=mol_docked.coordinates,
        atomicn1=mol_ref.atomicnums,
        atomicn2=mol_docked.atomicnums,
        minimize=minimize
    )
    return RMSD