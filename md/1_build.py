import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
import sys

pdbfile = '2lzm.pdb'

def prepare_protein(pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.

    CC-BY 4.0, TeachOpenCADD T019,  Molecular dynamics simulation
    (c) Pietro Gerletti, Mareike Leja, Jeffrey R Wagner, David Schaller, Andrea Volkamer, 2020.
    https://projects.volkamerlab.org/teachopencadd/talktorials/T019_md_simulation.html
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer


# Load the PDB structure
fixer = prepare_protein(pdbfile, ignore_missing_residues=False, ph=7.0)

# Solvate
modeller = app.Modeller(fixer.topology, fixer.positions)
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)

# Save topology and positions
with open('system.pdb', 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

# System Configuration
nonbondedMethod = app.PME
nonbondedCutoff = 1.0*unit.nanometers
ewaldErrorTolerance = 0.0005
constraints = app.HBonds
rigidWater = True
hydrogenMass = 1.5*unit.amu

# Create system
system = forcefield.createSystem(modeller.topology,
                                 nonbondedMethod=nonbondedMethod,
                                 nonbondedCutoff=nonbondedCutoff,
                                 constraints=constraints,
                                 rigidWater=rigidWater,
                                 ewaldErrorTolerance=ewaldErrorTolerance,
                                 hydrogenMass=hydrogenMass)

# Save system
with open('system.xml', 'w') as f:
    f.write(mm.openmm.XmlSerializer.serialize(system))

