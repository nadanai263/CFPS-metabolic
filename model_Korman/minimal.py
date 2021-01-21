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
	m = y[8]
	P = y[9]
	Pmat = y[10]

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
	Vmaxtx = params[12]
	d = params[13]
	KTX = params[14]
	Katp = params[15]
	kdeg = params[16]
	Vmaxtl = params[17]
	KTL = params[18]
	kmat = params[19]
	Vmax = params[20]

	derivs = [
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) - 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp),
	- 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) + 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) + 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp),
	+ 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	- 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	+ 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1) + 1.0 * 0.001*pep + 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) + 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp),
	- 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	- 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 1.0 * 0.001*pep,
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	+ 1.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) - 1.0 * kdeg*m,
	+ 1.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp) - 1.0 * kmat*P,
	+ 1.0 * kmat*P]
	return derivs

keysVar = ['atp','adp','nad','nadh','ppi','pppi','pep','pyr','m','P','Pmat']
valuesVar = [2.0,0.0,0.25,0.0,25.0,0.0,33.0,0.0,0.0,0.0,0.0]
dictVar = dict(zip(keysVar,valuesVar))

keysPar = ['VmaxPyk','Keq11','Kadp11','Kpep11','Katp11','Kpyr11','VmaxPpase','Keq22','Kpppi22','Kppi22','VmaxNoxe','Knadh23','Vmaxtx','d','KTX','Katp','kdeg','Vmaxtl','KTL','kmat','Vmax']
valuesPar = [1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1,2e-05,1e-06,1e-05,1.0,0.0001,2e-06,0.0001,0.0005,1.0]
dictPar = dict(zip(keysPar,valuesPar))

