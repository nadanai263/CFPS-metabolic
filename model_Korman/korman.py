# Python model ODEs from antimony file

import numpy as np

def model(y, t, params):

	atp = y[0]
	adp = y[1]
	nad = y[2]
	nadh = y[3]
	nadp = y[4]
	nadph = y[5]
	co2 = y[6]
	ppi = y[7]
	pppi = y[8]
	coa = y[9]
	glc = y[10]
	g6p = y[11]
	f6p = y[12]
	fbp = y[13]
	gap = y[14]
	dhap = y[15]
	bpg = y[16]
	p3g = y[17]
	p2g = y[18]
	pep = y[19]
	pyr = y[20]
	accoa = y[21]
	acaccoa = y[22]
	hmgcoa = y[23]
	mev = y[24]
	mvp = y[25]
	mpp = y[26]
	ipp = y[27]
	dmapp = y[28]
	gpp = y[29]
	lim = y[30]

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
	VmaxmGap = params[32]
	Keq7 = params[33]
	Kgap7 = params[34]
	Knadp7 = params[35]
	Kppi7 = params[36]
	Kbpg7 = params[37]
	Knadph7 = params[38]
	VmaxPgk = params[39]
	Keq8 = params[40]
	Kbpg8 = params[41]
	Kadp8 = params[42]
	Kp3g8 = params[43]
	Katp8 = params[44]
	VmaxPgm = params[45]
	Keq9 = params[46]
	Kp3g9 = params[47]
	Kp2g9 = params[48]
	VmaxEno = params[49]
	Keq10 = params[50]
	Kp2g10 = params[51]
	Kpep10 = params[52]
	VmaxPyk = params[53]
	Keq11 = params[54]
	Kadp11 = params[55]
	Kpep11 = params[56]
	Katp11 = params[57]
	Kpyr11 = params[58]
	VmaxPdh = params[59]
	Keq12 = params[60]
	Kcoa12 = params[61]
	Knad12 = params[62]
	Kpyr12 = params[63]
	Kaccoa12 = params[64]
	Kco212 = params[65]
	Knadh12 = params[66]
	VmaxThl = params[67]
	Keq13 = params[68]
	Kaccoa13 = params[69]
	Kacaccoa13 = params[70]
	Kcoa13 = params[71]
	VmaxHmgs = params[72]
	Keq14 = params[73]
	Kacaccoa14 = params[74]
	Kaccoa14 = params[75]
	Khmgcoa14 = params[76]
	Kcoa14 = params[77]
	VmaxHmgr = params[78]
	Keq15 = params[79]
	Khmgcoa15 = params[80]
	Knadph15 = params[81]
	Kmev15 = params[82]
	Kcoa15 = params[83]
	Knadp15 = params[84]
	VmaxMvk = params[85]
	Keq16 = params[86]
	Kmev16 = params[87]
	Katp16 = params[88]
	Kmvp16 = params[89]
	Kadp16 = params[90]
	VmaxPmvk = params[91]
	Keq17 = params[92]
	Kmvp17 = params[93]
	Katp17 = params[94]
	Kmpp17 = params[95]
	Kadp17 = params[96]
	VmaxMdc = params[97]
	Keq18 = params[98]
	Kmpp18 = params[99]
	Katp18 = params[100]
	Kadp18 = params[101]
	Kco218 = params[102]
	Kipp18 = params[103]
	Kppi18 = params[104]
	VmaxIdi = params[105]
	Keq19 = params[106]
	Kipp19 = params[107]
	Kdmapp19 = params[108]
	VmaxFS82F = params[109]
	Keq20 = params[110]
	Kdmapp20 = params[111]
	Kipp20 = params[112]
	Kgpp20 = params[113]
	Kpppi20 = params[114]
	VmaxLim = params[115]
	Keq21 = params[116]
	Kgpp21 = params[117]
	Klim21 = params[118]
	Kpppi21 = params[119]
	VmaxPpase = params[120]
	Keq22 = params[121]
	Kpppi22 = params[122]
	Kppi22 = params[123]
	VmaxNoxe = params[124]
	Knadh23 = params[125]

	derivs = [
	- 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) + 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 1.0 * VmaxMvk*(mev*atp-mvp*adp/Keq16)/(Kmev16*Katp16)*1/((1+mev/Kmev16+mvp/Kmvp16)*(1+atp/Katp16+adp/Kadp16)) - 1.0 * VmaxPmvk*(mvp*atp-mpp*adp/Keq17)/(Kmvp17*Katp17)*1/((1+mvp/Kmvp17+mpp/Kmpp17)*(1+atp/Katp17+adp/Kadp17)) - 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1),
	+ 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) + 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) + 1.0 * VmaxMvk*(mev*atp-mvp*adp/Keq16)/(Kmev16*Katp16)*1/((1+mev/Kmev16+mvp/Kmvp16)*(1+atp/Katp16+adp/Kadp16)) + 1.0 * VmaxPmvk*(mvp*atp-mpp*adp/Keq17)/(Kmvp17*Katp17)*1/((1+mvp/Kmvp17+mpp/Kmpp17)*(1+atp/Katp17+adp/Kadp17)) + 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1),
	- 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1) + 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	+ 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1) - 4.0 * VmaxNoxe*nadh/(Knadh23+nadh),
	- 2.0 * VmaxmGap*(gap*nadp*ppi-bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7)*1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1) + 2.0 * VmaxHmgr*(hmgcoa*nadph*nadph-mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15)*1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1),
	+ 2.0 * VmaxmGap*(gap*nadp*ppi-bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7)*1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1) - 2.0 * VmaxHmgr*(hmgcoa*nadph*nadph-mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15)*1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1),
	+ 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1) + 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1),
	- 3.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 3.0 * VmaxmGap*(gap*nadp*ppi-bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7)*1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1) + 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1) + 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	+ 1.0 * VmaxFS82F*(dmapp*ipp-gpp*pppi/Keq20)/(Kdmapp20*Kipp20)*1/((1+dmapp/Kdmapp20)*(1+ipp/Kipp20)+(1+gpp/Kgpp20)*(1+pppi/Kpppi20)-1) + 1.0 * VmaxLim*(gpp-lim*pppi/Keq21)/Kgpp21*1/(1+gpp/Kgpp21+(1+lim/Klim21)*(1+pppi/Kpppi21)-1) - 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1),
	- 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1) + 1.0 * VmaxThl*(accoa*accoa-coa*acaccoa/Keq13)/(Kaccoa13*Kaccoa13)*1/((1+accoa/Kaccoa13)*(1+accoa/Kaccoa13)+(1+acaccoa/Kacaccoa13)*(1+coa/Kcoa13)-1) + 1.0 * VmaxHmgs*(acaccoa*accoa-coa*hmgcoa/Keq14)/(Kacaccoa14*Kaccoa14)*1/((1+acaccoa/Kacaccoa14)*(1+accoa/Kaccoa14)+(1+hmgcoa/Khmgcoa14)*(1+coa/Kcoa14)-1) + 1.0 * VmaxHmgr*(hmgcoa*nadph*nadph-mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15)*1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1),
	- 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1),
	+ 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1),
	+ 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1),
	+ 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1),
	+ 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) + 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1) - 3.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 3.0 * VmaxmGap*(gap*nadp*ppi-bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7)*1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1),
	+ 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) - 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1),
	+ 3.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 3.0 * VmaxmGap*(gap*nadp*ppi-bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7)*1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1),
	+ 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1),
	+ 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1) - 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1),
	+ 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1),
	+ 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1),
	+ 3.0 * VmaxPdh*(coa*nad*pyr-accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12)*1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1) - 2.0 * VmaxThl*(accoa*accoa-coa*acaccoa/Keq13)/(Kaccoa13*Kaccoa13)*1/((1+accoa/Kaccoa13)*(1+accoa/Kaccoa13)+(1+acaccoa/Kacaccoa13)*(1+coa/Kcoa13)-1) - 1.0 * VmaxHmgs*(acaccoa*accoa-coa*hmgcoa/Keq14)/(Kacaccoa14*Kaccoa14)*1/((1+acaccoa/Kacaccoa14)*(1+accoa/Kaccoa14)+(1+hmgcoa/Khmgcoa14)*(1+coa/Kcoa14)-1),
	+ 1.0 * VmaxThl*(accoa*accoa-coa*acaccoa/Keq13)/(Kaccoa13*Kaccoa13)*1/((1+accoa/Kaccoa13)*(1+accoa/Kaccoa13)+(1+acaccoa/Kacaccoa13)*(1+coa/Kcoa13)-1) - 1.0 * VmaxHmgs*(acaccoa*accoa-coa*hmgcoa/Keq14)/(Kacaccoa14*Kaccoa14)*1/((1+acaccoa/Kacaccoa14)*(1+accoa/Kaccoa14)+(1+hmgcoa/Khmgcoa14)*(1+coa/Kcoa14)-1),
	+ 1.0 * VmaxHmgs*(acaccoa*accoa-coa*hmgcoa/Keq14)/(Kacaccoa14*Kaccoa14)*1/((1+acaccoa/Kacaccoa14)*(1+accoa/Kaccoa14)+(1+hmgcoa/Khmgcoa14)*(1+coa/Kcoa14)-1) - 1.0 * VmaxHmgr*(hmgcoa*nadph*nadph-mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15)*1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1),
	+ 1.0 * VmaxHmgr*(hmgcoa*nadph*nadph-mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15)*1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1) - 1.0 * VmaxMvk*(mev*atp-mvp*adp/Keq16)/(Kmev16*Katp16)*1/((1+mev/Kmev16+mvp/Kmvp16)*(1+atp/Katp16+adp/Kadp16)),
	+ 1.0 * VmaxMvk*(mev*atp-mvp*adp/Keq16)/(Kmev16*Katp16)*1/((1+mev/Kmev16+mvp/Kmvp16)*(1+atp/Katp16+adp/Kadp16)) - 1.0 * VmaxPmvk*(mvp*atp-mpp*adp/Keq17)/(Kmvp17*Katp17)*1/((1+mvp/Kmvp17+mpp/Kmpp17)*(1+atp/Katp17+adp/Kadp17)),
	+ 1.0 * VmaxPmvk*(mvp*atp-mpp*adp/Keq17)/(Kmvp17*Katp17)*1/((1+mvp/Kmvp17+mpp/Kmpp17)*(1+atp/Katp17+adp/Kadp17)) - 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1),
	+ 1.0 * VmaxMdc*(mpp*atp-adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18)*1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1) - 1.0 * VmaxIdi*(ipp-dmapp/Keq19)/Kipp19*1/(1+ipp/Kipp19+1+dmapp/Kdmapp19-1) - 1.0 * VmaxFS82F*(dmapp*ipp-gpp*pppi/Keq20)/(Kdmapp20*Kipp20)*1/((1+dmapp/Kdmapp20)*(1+ipp/Kipp20)+(1+gpp/Kgpp20)*(1+pppi/Kpppi20)-1),
	+ 1.0 * VmaxIdi*(ipp-dmapp/Keq19)/Kipp19*1/(1+ipp/Kipp19+1+dmapp/Kdmapp19-1) - 1.0 * VmaxFS82F*(dmapp*ipp-gpp*pppi/Keq20)/(Kdmapp20*Kipp20)*1/((1+dmapp/Kdmapp20)*(1+ipp/Kipp20)+(1+gpp/Kgpp20)*(1+pppi/Kpppi20)-1),
	+ 1.0 * VmaxFS82F*(dmapp*ipp-gpp*pppi/Keq20)/(Kdmapp20*Kipp20)*1/((1+dmapp/Kdmapp20)*(1+ipp/Kipp20)+(1+gpp/Kgpp20)*(1+pppi/Kpppi20)-1) - 1.0 * VmaxLim*(gpp-lim*pppi/Keq21)/Kgpp21*1/(1+gpp/Kgpp21+(1+lim/Klim21)*(1+pppi/Kpppi21)-1),
	+ 1.0 * VmaxLim*(gpp-lim*pppi/Keq21)/Kgpp21*1/(1+gpp/Kgpp21+(1+lim/Klim21)*(1+pppi/Kpppi21)-1)]
	return derivs

keysVar = ['atp','adp','nad','nadh','nadp','nadph','co2','ppi','pppi','coa','glc','g6p','f6p','fbp','gap','dhap','bpg','p3g','p2g','pep','pyr','accoa','acaccoa','hmgcoa','mev','mvp','mpp','ipp','dmapp','gpp','lim']
valuesVar = [4.0,0.0,0.25,0.0,1.5,0.0,0.0,25.0,0.0,1.5,500.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
dictVar = dict(zip(keysVar,valuesVar))

keysPar = ['VmaxHex','Keq1','Katp1','Kglc1','Kadp1','Kg6p1','VmaxPgi','Keq2','Kg6p2','Kf6p2','VmaxPfk','Keq3','Katp3','Kf6p3','Kadp3','Kfbp3','VmaxFba','Keq4','Kfbp4','Kdhap4','Kgap4','VmaxTpi','Keq5','Kdhap5','Kgap5','VmaxGap','Keq6','Kgap6','Knad6','Kppi6','Kbpg6','Knadh6','VmaxmGap','Keq7','Kgap7','Knadp7','Kppi7','Kbpg7','Knadph7','VmaxPgk','Keq8','Kbpg8','Kadp8','Kp3g8','Katp8','VmaxPgm','Keq9','Kp3g9','Kp2g9','VmaxEno','Keq10','Kp2g10','Kpep10','VmaxPyk','Keq11','Kadp11','Kpep11','Katp11','Kpyr11','VmaxPdh','Keq12','Kcoa12','Knad12','Kpyr12','Kaccoa12','Kco212','Knadh12','VmaxThl','Keq13','Kaccoa13','Kacaccoa13','Kcoa13','VmaxHmgs','Keq14','Kacaccoa14','Kaccoa14','Khmgcoa14','Kcoa14','VmaxHmgr','Keq15','Khmgcoa15','Knadph15','Kmev15','Kcoa15','Knadp15','VmaxMvk','Keq16','Kmev16','Katp16','Kmvp16','Kadp16','VmaxPmvk','Keq17','Kmvp17','Katp17','Kmpp17','Kadp17','VmaxMdc','Keq18','Kmpp18','Katp18','Kadp18','Kco218','Kipp18','Kppi18','VmaxIdi','Keq19','Kipp19','Kdmapp19','VmaxFS82F','Keq20','Kdmapp20','Kipp20','Kgpp20','Kpppi20','VmaxLim','Keq21','Kgpp21','Klim21','Kpppi21','VmaxPpase','Keq22','Kpppi22','Kppi22','VmaxNoxe','Knadh23']
valuesPar = [0.05,1000.0,0.1,0.15,0.1,0.1,1.0,0.5,0.3,0.15,0.5,300.0,0.1,0.03,0.1,0.1,1.0,0.2,0.015,0.1,0.1,0.25,0.04,1.0,1.0,0.5,0.07,1.0,0.05,0.5,2.0,0.1,10.0,0.07,1.0,0.05,0.5,2.0,0.1,2.0,3200.0,0.1,0.1,2.2,3.0,1.0,0.19,0.5,0.1,1.0,6.7,0.1,1.0,1.0,5000.0,0.3,0.3,1.0,1.0,10.0,5000.0,0.1,0.1,0.5,0.1,0.1,0.1,1.0,0.05,0.4,0.1,0.1,1.0,50.0,0.01,0.4,0.5,0.5,1.0,50.0,0.02,0.1,1.0,0.3,0.3,5.0,10.0,0.07,0.5,1.5,0.1,5.0,2.0,0.008,0.14,0.01,0.4,5.0,10.0,0.1,0.06,0.3,5.0,0.5,0.5,1.0,0.77,0.004,0.05,1.0,200.0,0.005,0.005,4.0,5.0,5.0,5000.0,0.002,2.0,0.5,20.0,10.0,0.1,0.1,1.0,0.1]
dictPar = dict(zip(keysPar,valuesPar))

