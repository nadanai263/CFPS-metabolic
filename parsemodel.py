# Load model from Antimony or SBML file
# Extract ODEs and generate Python/Julia scripts for simulation
# Save a SBML file if required

import sbmlIO as sio
import tellurium as te
import os

#here = "/app" # Docker image mount point. 
here = os.getcwd() # Use this if running locally.

PATH_TO_MODELS = here+"/model_Korman/"
PATH_TO_OUTPUT = here+"/modeloutput/" # Path for output data, plots
MODELNAME = "autocatalytic"
MODELPY = "autocatalytic.py"
MODELJU = "autocatalytic.jl"
MODELSBML = "autocatalytic_SBML.xml"
FILENAME = "autocatalytic.csv" # Filename for output data

# Generate ODEs 
r = te.loada(PATH_TO_MODELS+MODELNAME) # from Antimony file
#r = te.loadSBMLModel(PATH_TO_MODELS+MODELexSBML) # or from SBML
odes = sio.getODEsFromModel(r)
speciesIds, speciesValues, parameterIds, parameterValues, derivatives = sio.parseODEs(r,odes)

print(speciesIds)
print(speciesValues)
print(parameterIds)
print(parameterValues)
sio.writePython(speciesIds,speciesValues,parameterIds,parameterValues,derivatives,PATH_TO_MODELS,MODELPY)
sio.writeJulia(speciesIds,speciesValues,parameterIds,parameterValues,derivatives,PATH_TO_MODELS,MODELJU)
print(odes)

# Save model as SBML
r.exportToSBML(PATH_TO_MODELS+MODELSBML)

