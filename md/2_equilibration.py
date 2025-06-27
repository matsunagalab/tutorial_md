import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
import sys

# Load the topology and positions
with open('system.pdb', 'r') as f:
    pdb = app.PDBFile(f)

# Load the system
with open('system.xml', 'r') as f:
    system = mm.openmm.XmlSerializer.deserialize(f.read())

# Integrator options
constraintTolerance = 0.000001
dt = 0.004*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1.0/unit.picosecond
pressure = 1.0*unit.atmospheres
barostatInterval = 25

# Simulation options
steps = 25000

# Linux with GPU
#platform = mm.openmm.Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}

# Mac
platform = mm.openmm.Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'single'}

dcdReporter = app.DCDReporter('2_equilibration.dcd', 100)
dataReporter = app.StateDataReporter(sys.stdout, 100, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
#checkpointReporter = CheckpointReporter('run1_checkpoint.chk', 10000)

# Prepare simulation object
print('Building system...')
system.addForce(mm.openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = mm.openmm.LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)

# Positional restraints
restraint = mm.openmm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
#restraint.addGlobalParameter('k', 10.0*unit.kilojoules_per_mole/(unit.nanometer*unit.nanometer))
restraint.addPerParticleParameter('k')
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')
for atom in pdb.topology.atoms():
    if atom.name == 'CA':
        restraint.addParticle(atom.index, [100.0*unit.kilojoules_per_mole/(unit.nanometer*unit.nanometer), 
              pdb.positions[atom.index][0], pdb.positions[atom.index][1], pdb.positions[atom.index][2]] )

simulation = app.Simulation(pdb.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# Minimize
print('Performing energy minimization...')
simulation.minimizeEnergy()

# Equilibrate
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.step(steps)

#system.removeForce(system.getNumForces() - 1)
simulation.saveState("2_equilibration.xml")

