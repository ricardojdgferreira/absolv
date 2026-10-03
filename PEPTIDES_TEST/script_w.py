# general imports
import os, sys, re
import pandas as pd

# import absolv functions
import absolv.config
import absolv.runner
from absolv.density import Density

# import openmm functions
import openmm.unit as unit
import femto.md.config
import openmm.unit
from openmm.app import *

# import openff toolkit
import openff.toolkit

# helper functions
def mol_weight_smiles(smiles):
	mol = Chem.MolFromSmiles(smiles)
	return Descriptors.MolWt(mol)

def mol_weight_sdf(sdf):
	mol = Chem.MolFromMolFile(sdf)
	return Descriptors.MolWt(mol)

def component_density(smiles):
	if smiles == "O": # water
		return 0.997 				# g/cm3 @ 298 K
	return float(Density(smiles))   # your density predictor

# parameters
name = "peptide37_w"
solutes={"peptide_37.sdf": 1}
solvent_a=None
solvent_b={"O": 600}

#get density for solvents
if solvent_a == None:
	box_target_density_a=None
else:
	smiles_a = next(iter(solvent_a))
	box_target_density_a = component_density(smiles_a)
	
smiles_b = next(iter(solvent_b))
box_target_density_b = component_density(smiles_b)

if box_target_density_a is None:
	print(f"Target density for solvent_a is {box_target_density_a} (vacuum)")
else:
	print(f"Target density for solvent_a is {box_target_density_a} g/ml")
print(f"Target density for solvent_b is {box_target_density_b} g/ml")

# define system
system=absolv.config.System(solutes=solutes, solvent_a=solvent_a, solvent_b=solvent_b)

# set temperature and pressure
temperature = 298.15 * openmm.unit.kelvin
pressure = 1.0 * openmm.unit.atmosphere

# define non-eq protocol
integrator = femto.md.config.LangevinIntegrator(timestep = 2.0 * openmm.unit.femtoseconds)

alchemical_protocol_a=absolv.config.NonEquilibriumProtocol(
	# the protocol to use for the production run at each end state
	production_protocol=absolv.config.SimulationProtocol(
		integrator=integrator, n_steps=12500 * 160 # 2 ns (6250 * 160 = 1.000.000 steps)
	),
	# the interval with which to store frames that NEQ switching simulations will be launched from
	production_report_interval=6250,
	# define how the NEQ switching will be performed
	switching_protocol=absolv.config.SwitchingProtocol(
		# Annihilate the electrostatic interactions over the first 12 ps
		n_electrostatic_steps=60,
		n_steps_per_electrostatic_step=100,
		# followed by decoupling the vdW interactions over the next 38 ps (both are zero if vaccum to be used)
		n_steric_steps=0,
		n_steps_per_steric_step=0
		)
	)

alchemical_protocol_b=absolv.config.NonEquilibriumProtocol(
    production_protocol=absolv.config.SimulationProtocol(
        integrator=integrator, n_steps=25000 * 160,
	production_report_interval=12500,
	),
	switching_protocol=absolv.config.SwitchingProtocol(
		# Annihilate the electrostatic interactions over the first 12 ps
		n_electrostatic_steps=60,
		n_steps_per_electrostatic_step=100,
		# followed by decoupling the vdW interactions over the next 38 ps
		n_steric_steps=190,
		n_steps_per_steric_step=100,
	    )
    )

# get components
config = absolv.config.Config(temperature=temperature,pressure=pressure,alchemical_protocol_a=alchemical_protocol_a,alchemical_protocol_b=alchemical_protocol_b,)

# define forcefield and setup runs
force_field = openff.toolkit.ForceField("openff-2.2.0.offxml")
if box_target_density_a is None:
    prepared_system_a, prepared_system_b = absolv.runner.setup(system, config, force_field,
                                                  box_target_density_b=(box_target_density_b * openmm.unit.grams / openmm.unit.milliliters),)
else:
   prepared_system_a, prepared_system_b = absolv.runner.setup(system, config, force_field,
                                                  box_target_density_a=(box_target_density_a * openmm.unit.grams / openmm.unit.milliliters),
                                                  box_target_density_b=(box_target_density_b * openmm.unit.grams / openmm.unit.milliliters),)

# save coordinates
path_to_save = os.path.join(os.getcwd(), name)
os.system("mkdir " + path_to_save)

with open(path_to_save + '/system_a.pdb', 'w') as output_file:
	PDBFile.writeFile(prepared_system_a.topology, prepared_system_a.coords, output_file)
with open(path_to_save + '/system_b.pdb', 'w') as output_file:
	PDBFile.writeFile(prepared_system_b.topology, prepared_system_b.coords, output_file)

# run simulations
result = absolv.runner.run_neq(config, prepared_system_a, prepared_system_b, "CUDA", path_to_save)

# save results
with open(path_to_save + '/results.txt', 'w') as output_file:
	for i in result:
		output_file.write(str(i) + '\n')
