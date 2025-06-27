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
#steps = 250000000
steps = 250000

# Linux with GPU
#platform = mm.openmm.Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}

# Mac
platform = mm.openmm.Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'single'}

dcdReporter = app.DCDReporter('3_production.dcd', 250)
dataReporter = app.StateDataReporter(sys.stdout, 250, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
#checkpointReporter = CheckpointReporter('run1_checkpoint.chk', 10000)

# Prepare simulation object
print('Building system...')
system.addForce(mm.openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = mm.openmm.LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = app.Simulation(pdb.topology, system, integrator, platform, platformProperties)
#simulation.context.setPositions(pdb.positions)

# Load the state
simulation.loadState('2_equilibration.xml')

# Production
print('Production...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)

simulation.saveState("3_production.xml")

