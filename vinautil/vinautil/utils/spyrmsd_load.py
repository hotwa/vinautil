from spyrmsd.molecule import Molecule
from spyrmsd import io, rmsd
from spyrmsd.exceptions import NonIsomorphicGraphs
from numpy import nan as np_nan
from spyrmsd.optional.obabel import to_molecule
from openbabel import pybel
from typing import List

def symmrmsd_mol2_list(mol2_ref: str, mol2_docked: List[str], minimize=False) -> float:
    ref = to_molecule(pybel.readstring('mol2',mol2_ref))
    mol_lst = [to_molecule(pybel.readstring('mol2',i)) for i in mol2_docked]
    mol = mol_lst[0]
    ref.strip()  # remove Hydrogen atom
    mol.strip()  # remove Hydrogen atom
    coords_ref = ref.coordinates
    anum_ref = ref.atomicnums
    adj_ref = ref.adjacency_matrix
    coords = [mol.coordinates for mol in mol_lst]
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
        return RMSD
    except NonIsomorphicGraphs:
        return np_nan
    
def rmsd_mol2_list(mol2_ref: str, mol2_docked: List[str], minimize=False) -> List[float]:
    ref = to_molecule(pybel.readstring('mol2',mol2_ref))
    ref.strip()  # remove Hydrogen atom
    coords_ref = ref.coordinates
    mol_ref_atomnum = ref.atomicnums

    RMSDs = []
    for mol2 in mol2_docked:
        try:
            mol = to_molecule(pybel.readstring('mol2',mol2))
            mol.strip()  # remove Hydrogen atom
            RMSD = rmsd.rmsd(
                coords1=coords_ref,
                coords2=mol.coordinates,
                atomicn1=mol_ref_atomnum,
                atomicn2=mol.atomicnums,
                minimize=minimize
            )
            RMSDs.append(RMSD)
        except:
            RMSDs.append(np_nan)
    return RMSDs

def symmrmsd_mol2(mol2_ref:str, mol2_docked:str, minimize=False) -> float:
    ref = to_molecule(pybel.readstring(mol2_ref, 'mol2'))
    mol = to_molecule(pybel.readstring(mol2_docked, 'mol2'))
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

def rmsd_mol2(mol2_ref: str, mol2_docked: str, minimize=False) -> float:
    mol_ref = to_molecule(pybel.readstring(mol2_ref, 'mol2'))
    mol_docked = to_molecule(pybel.readstring(mol2_docked, 'mol2'))
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