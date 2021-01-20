# Julia model ODEs from antimony file

function model!(du, u, params, t)

	atp,adp,nad,nadh,ppi,pppi,pep,pyr = u

	VmaxPyk,Keq11,Kadp11,Kpep11,Katp11,Kpyr11,VmaxPpase,Keq22,Kpppi22,Kppi22,VmaxNoxe,Knadh23,Vmax = params

	du[1] = + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) - 1.0 * Vmax*atp
	du[2] = - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1) + 1.0 * Vmax*atp
	du[3] = + 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[4] = - 4.0 * VmaxNoxe*nadh/(Knadh23+nadh)
	du[5] = + 4.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1) + 1.0 * Vmax*atp
	du[6] = - 2.0 * VmaxPpase*(pppi-ppi*ppi/Keq22)/Kpppi22*1/(1+pppi/Kpppi22+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1)
	du[7] = - 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
	du[8] = + 3.0 * VmaxPyk*(pep*adp-pyr*atp/Keq11)/(Kadp11*Kpep11)*1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)
end

keysVar = ["atp","adp","nad","nadh","ppi","pppi","pep","pyr"]
valuesVar = [4.0,0.0,0.25,0.0,25.0,0.0,0.0,0.0]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["VmaxPyk","Keq11","Kadp11","Kpep11","Katp11","Kpyr11","VmaxPpase","Keq22","Kpppi22","Kppi22","VmaxNoxe","Knadh23","Vmax"]
valuesPar = [1.0,5000.0,0.3,0.3,1.0,1.0,20.0,10.0,0.1,0.1,1.0,0.1,1.0]
dictPar = Dict(keysPar .=> valuesPar)

