# Python model ODEs from antimony file

import numpy as np

def model(y, t, params):

	atp = y[0]
	adp = y[1]
	nad = y[2]
	nadh = y[3]
	ppi = y[4]
	pppi = y[5]
	glc = y[6]
	g6p = y[7]
	f6p = y[8]
	fbp = y[9]
	gap = y[10]
	dhap = y[11]
	bpg = y[12]
	p3g = y[13]
	p2g = y[14]
	pep = y[15]
	pyr = y[16]

	VmaxHex = params[0]
	Keq1 = params[1]
	Katp1 = params[2]
	Kglc1 = params[3]
	Kadp1 = params[4]
	Kg6p1 = params[5]
	VmaxPgi = params[6]
	Keq2 = params[7]
	Kg6p2 = params[8]
	Kf6p2 = params[9]
	VmaxPfk = params[10]
	Keq3 = params[11]
	Katp3 = params[12]
	Kf6p3 = params[13]
	Kadp3 = params[14]
	Kfbp3 = params[15]
	VmaxFba = params[16]
	Keq4 = params[17]
	Kfbp4 = params[18]
	Kdhap4 = params[19]
	Kgap4 = params[20]
	VmaxTpi = params[21]
	Keq5 = params[22]
	Kdhap5 = params[23]
	Kgap5 = params[24]
	VmaxGap = params[25]
	Keq6 = params[26]
	Kgap6 = params[27]
	Knad6 = params[28]
	Kppi6 = params[29]
	Kbpg6 = params[30]
	Knadh6 = params[31]
	VmaxPgk = params[32]
	Keq8 = params[33]
	Kbpg8 = params[34]
	Kadp8 = params[35]
	Kp3g8 = params[36]
	Katp8 = params[37]
	VmaxPgm = params[38]
	Keq9 = params[39]
	Kp3g9 = params[40]
	Kp2g9 = params[41]
	VmaxEno = params[42]
	Keq10 = params[43]
	Kp2g10 = params[44]
	Kpep10 = params[45]
	VmaxPyk = params[46]
	Keq11 = params[47]
	Kadp11 = params[48]
	Kpep11 = params[49]
	Katp11 = params[50]
	Kpyr11 = params[51]
	VmaxPpase = params[52]
	Keq22 = params[53]
	Kpppi22 = params[54]
	Kppi22 = params[55]
	VmaxNoxe = params[56]
	Knadh23 = params[57]

	derivs = [
	- 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) + 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	+ 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) + 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	- 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	+ 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	- 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	- 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	- 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1),
	+ 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1),
	+ 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1),
	+ 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1),
	+ 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) + 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1) - 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1),
	+ 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) - 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1),
	+ 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1),
	+ 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1),
	+ 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1) - 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1),
	+ 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)]
	return derivs

keysVar = ['atp','adp','nad','nadh','ppi','pppi','glc','g6p','f6p','fbp','gap','dhap','bpg','p3g','p2g','pep','pyr']
valuesVar = [4.0,0.0,0.25,0.0,25.0,0.0,500.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
dictVar = dict(zip(keysVar,valuesVar))

keysPar = ['VmaxHex','Keq1','Katp1','Kglc1','Kadp1','Kg6p1','VmaxPgi','Keq2','Kg6p2','Kf6p2','VmaxPfk','Keq3','Katp3','Kf6p3','Kadp3','Kfbp3','VmaxFba','Keq4','Kfbp4','Kdhap4','Kgap4','VmaxTpi','Keq5','Kdhap5','Kgap5','VmaxGap','Keq6','Kgap6','Knad6','Kppi6','Kbpg6','Knadh6','VmaxPgk','Keq8','Kbpg8','Kadp8','Kp3g8','Katp8','VmaxPgm','Keq9','Kp3g9','Kp2g9','VmaxEno','Keq10','Kp2g10','Kpep10','VmaxPyk','Keq11','Kadp11','Kpep11','Katp11','Kpyr11','VmaxPpase','Keq22','Kpppi22','Kppi22','VmaxNoxe','Knadh23']
valuesPar = [0.05,1000.0,0.1,0.15,0.1,0.1,1.0,0.5,0.3,0.15,0.5,300.0,0.1,0.03,0.1,0.1,1.0,0.2,0.015,0.1,0.1,0.25,0.04,1.0,1.0,0.5,0.07,1.0,0.05,0.5,2.0,0.1,2.0,3200.0,0.1,0.1,2.2,3.0,1.0,0.19,0.5,0.1,1.0,6.7,0.1,1.0,1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1]
dictPar = dict(zip(keysPar,valuesPar))

