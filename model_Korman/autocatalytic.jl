# Julia model ODEs from antimony file

function model!(du, u, params, t)

	A,B,M,E,R,W = u

	kcat1,Ka,Kb,Km,dG0,T,kcat2,Kr,Kw,Ke,dG1,kdeg = params

	du[1] = - 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A))))
	du[2] = + 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A)))) - 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M)))))
	du[3] = + 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A)))) - 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M)))))
	du[4] = + 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) - 1.0 * (kdeg*E)
	du[5] = + 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) - 1.0 * (kdeg*R)
	du[6] = + 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) + 1.0 * (kdeg*R) + 1.0 * (kdeg*E)
end

keysVar = ["A","B","M","E","R","W"]
valuesVar = [50.0,0.1,0.1,0.05,8.3145,0.1]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["kcat1","Ka","Kb","Km","dG0","T","kcat2","Kr","Kw","Ke","dG1","kdeg"]
valuesPar = [10.0,0.1,0.1,0.1,-1000.0,298.0,10.0,0.1,0.1,0.1,-1000.0,0.1]
dictPar = Dict(keysPar .=> valuesPar)

