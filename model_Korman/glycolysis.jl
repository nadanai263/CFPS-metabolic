# Julia model ODEs from antimony file

function model!(du, u, params, t)

	atp,adp,nad,nadh,ppi,pppi,glc,g6p,f6p,fbp,gap,dhap,bpg,p3g,p2g,pep,pyr = u

	VmaxHex,Keq1,Katp1,Kglc1,Kadp1,Kg6p1,VmaxPgi,Keq2,Kg6p2,Kf6p2,VmaxPfk,Keq3,Katp3,Kf6p3,Kadp3,Kfbp3,VmaxFba,Keq4,Kfbp4,Kdhap4,Kgap4,VmaxTpi,Keq5,Kdhap5,Kgap5,VmaxGap,Keq6,Kgap6,Knad6,Kppi6,Kbpg6,Knadh6,VmaxPgk,Keq8,Kbpg8,Kadp8,Kp3g8,Katp8,VmaxPgm,Keq9,Kp3g9,Kp2g9,VmaxEno,Keq10,Kp2g10,Kpep10,VmaxPyk,Keq11,Kadp11,Kpep11,Katp11,Kpyr11,VmaxPpase,Keq22,Kpppi22,Kppi22,VmaxNoxe,Knadh23 = params

	du[1] = - 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) + 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
	du[2] = + 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) + 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
	du[3] = - 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[4] = + 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[5] = - 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) + 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1)
	du[6] = - 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1)
	du[7] = - 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1)
	du[8] = + 1.5 * VmaxHex*(atp*glc-g6p*adp/Keq1)/(Katp1*Kglc1)*1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1) - 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1)
	du[9] = + 1.5 * VmaxPgi*(g6p-f6p/Keq2)/Kg6p2*1/(1+g6p/Kg6p2+1+f6p/Kf6p2-1) - 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1)
	du[10] = + 1.5 * VmaxPfk*(f6p*atp-fbp*adp/Keq3)/(Katp3*Kf6p3)*1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1) - 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1)
	du[11] = + 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) + 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1) - 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1)
	du[12] = + 1.5 * VmaxFba*(fbp-gap*dhap/Keq4)/Kfbp4*1/(1+fbp/Kfbp4+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1) - 1.5 * VmaxTpi*(dhap-gap/Keq5)/Kdhap5*1/(1+dhap/Kdhap5+1+gap/Kgap5-1)
	du[13] = + 1.0 * VmaxGap*(gap*nad*ppi-bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6)*1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1) - 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1)
	du[14] = + 3.0 * VmaxPgk*(bpg*adp-p3g*atp/Keq8)/(Kbpg8*Kadp8)*1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1) - 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1)
	du[15] = + 3.0 * VmaxPgm*(p3g-p2g/Keq9)/Kp3g9*1/(1+p3g/Kp3g9+1+p2g/Kp2g9-1) - 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1)
	du[16] = + 3.0 * VmaxEno*(p2g-pep/Keq10)/Kp2g10*1/(1+p2g/Kp2g10+1+pep/Kpep10-1) - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
	du[17] = + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
end

keysVar = ["atp","adp","nad","nadh","ppi","pppi","glc","g6p","f6p","fbp","gap","dhap","bpg","p3g","p2g","pep","pyr"]
valuesVar = [4.0,0.0,0.25,0.0,25.0,0.0,500.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["VmaxHex","Keq1","Katp1","Kglc1","Kadp1","Kg6p1","VmaxPgi","Keq2","Kg6p2","Kf6p2","VmaxPfk","Keq3","Katp3","Kf6p3","Kadp3","Kfbp3","VmaxFba","Keq4","Kfbp4","Kdhap4","Kgap4","VmaxTpi","Keq5","Kdhap5","Kgap5","VmaxGap","Keq6","Kgap6","Knad6","Kppi6","Kbpg6","Knadh6","VmaxPgk","Keq8","Kbpg8","Kadp8","Kp3g8","Katp8","VmaxPgm","Keq9","Kp3g9","Kp2g9","VmaxEno","Keq10","Kp2g10","Kpep10","VmaxPyk","Keq11","Kadp11","Kpep11","Katp11","Kpyr11","VmaxPpase","Keq22","Kpppi22","Kppi22","VmaxNoxe","Knadh23"]
valuesPar = [0.05,1000.0,0.1,0.15,0.1,0.1,1.0,0.5,0.3,0.15,0.5,300.0,0.1,0.03,0.1,0.1,1.0,0.2,0.015,0.1,0.1,0.25,0.04,1.0,1.0,0.5,0.07,1.0,0.05,0.5,2.0,0.1,2.0,3200.0,0.1,0.1,2.2,3.0,1.0,0.19,0.5,0.1,1.0,6.7,0.1,1.0,1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1]
dictPar = Dict(keysPar .=> valuesPar)

