# Julia model ODEs from antimony file

function model!(du, u, params, t)

	atp,adp,nad,nadh,ppi,pppi,pep,pyr,m,P,Pmat = u

	VmaxPyk,Keq11,Kadp11,Kpep11,Katp11,Kpyr11,VmaxPpase,Keq22,Kpppi22,Kppi22,VmaxNoxe,Knadh23,Vmaxtx,d,KTX,Katp,kdeg,Vmaxtl,KTL,kmat,Vmax = params

	du[1] = + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) - 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp)
	du[2] = - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) + 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) + 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp)
	du[3] = + 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[4] = - 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[5] = + 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1) + 1.0 * 0.001*pep + 170.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) + 904.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp)
	du[6] = - 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1)
	du[7] = - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 1.0 * 0.001*pep
	du[8] = + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
	du[9] = + 1.0 * Vmaxtx*d/(d+KTX)*atp/(atp+Katp) - 1.0 * kdeg*m
	du[10] = + 1.0 * Vmaxtl*m/(m+KTL)*atp/(atp+Katp) - 1.0 * kmat*P
	du[11] = + 1.0 * kmat*P
end

keysVar = ["atp","adp","nad","nadh","ppi","pppi","pep","pyr","m","P","Pmat"]
valuesVar = [2.0,0.0,0.25,0.0,25.0,0.0,33.0,0.0,0.0,0.0,0.0]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["VmaxPyk","Keq11","Kadp11","Kpep11","Katp11","Kpyr11","VmaxPpase","Keq22","Kpppi22","Kppi22","VmaxNoxe","Knadh23","Vmaxtx","d","KTX","Katp","kdeg","Vmaxtl","KTL","kmat","Vmax"]
valuesPar = [1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1,2e-05,1e-06,1e-05,1.0,0.0001,2e-06,0.0001,0.0005,1.0]
dictPar = Dict(keysPar .=> valuesPar)

