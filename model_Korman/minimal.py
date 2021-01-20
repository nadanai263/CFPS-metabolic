# Python model ODEs from antimony file

import numpy as np

def model(y, t, params):

	atp = y[0]
	adp = y[1]
	nad = y[2]
	nadh = y[3]
	ppi = y[4]
	pppi = y[5]
	pep = y[6]
	pyr = y[7]

	VmaxPyk = params[0]
	Keq11 = params[1]
	Kadp11 = params[2]
	Kpep11 = params[3]
	Katp11 = params[4]
	Kpyr11 = params[5]
	VmaxPpase = params[6]
	Keq22 = params[7]
	Kpppi22 = params[8]
	Kppi22 = params[9]
	VmaxNoxe = params[10]
	Knadh23 = params[11]
	Vmax = params[12]

	derivs = [
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 1.0 * Vmax*atp,
	- 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) + 1.0 * Vmax*atp,
	+ 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	- 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	+ 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1) + 1.0 * Vmax*atp,
	- 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	- 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)]
	return derivs

keysVar = ['atp','adp','nad','nadh','ppi','pppi','pep','pyr']
valuesVar = [4.0,0.0,0.25,0.0,25.0,0.0,0.0,0.0]
dictVar = dict(zip(keysVar,valuesVar))

keysPar = ['VmaxPyk','Keq11','Kadp11','Kpep11','Katp11','Kpyr11','VmaxPpase','Keq22','Kpppi22','Kppi22','VmaxNoxe','Knadh23','Vmax']
valuesPar = [1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1,1.0]
dictPar = dict(zip(keysPar,valuesPar))

