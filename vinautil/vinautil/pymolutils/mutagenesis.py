# 定义残基突变函数
from pymol import cmd
from pathlib import Path
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
    cmd.reinitialize('everything')  # ! clean up
    cmd.load(filename.as_posix())
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
    cmd.reinitialize('everything')  # Reinitialize PyMOL to default settings.
    return _out_file