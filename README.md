## Install packages

```
$ git clone https://github.com/matsunagalab/tutorial_md.git
$ cd tutorial_md/
```

```
$ brew install uv
$ uv venv --python=python3.11
$ source .venv/bin/activate
$ uv pip install openmm mdtraj matplotlib scikit-learn numpy jupyterlab
$ uv pip install "git+https://github.com/openmm/pdbfixer.git@v1.11"
```

## Check Lysozyme structure on PDB

https://www.rcsb.org/structure/2LZM

```
Download the pdb file in md/ folder.
```

## Run MD with OpenMM

```
$ cd md/
$ python 1_build.py
$ python 2_equilibration.py
$ pytohn 3_production.py
```

```
# extract protein as protein.pdb and protein.dcd
$ python3 -c "import mdtraj as md; md.load('system.pdb', atom_indices=md.load('system.pdb').topology.select('protein')).save('protein.pdb')"
$ python3 -c "import mdtraj as md; md.load('3_production.dcd', top='system.pdb', atom_indices=md.load('system.pdb').topology.select('protein')).save('protein.dcd')"
# extract Ca atoms as ca.pdb and ca.dcd
$ python3 -c "import mdtraj as md; md.load('system.pdb', atom_indices=md.load('system.pdb').topology.select('name CA')).save('ca.pdb')"
$ python3 -c "import mdtraj as md; md.load('3_production.dcd', top='system.pdb', atom_indices=md.load('system.pdb').topology.select('name CA')).save('ca.dcd')"
```

## Analayze MD trajectory with MDTraj

```
$ cd analyze/
$ vscode ./ # or cursor ./
# open notebooks in VSCode or Cursor
```
