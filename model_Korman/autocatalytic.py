# Python model ODEs from antimony file

import numpy as np
from numpy import exp as exp
from numpy import log as log

def model(y, t, params):

	A = y[0]
	B = y[1]
	M = y[2]
	E = y[3]
	R = y[4]
	W = y[5]

	kcat1 = params[0]
	Ka = params[1]
	Kb = params[2]
	Km = params[3]
	dG0 = params[4]
	T = params[5]
	kcat2 = params[6]
	Kr = params[7]
	Kw = params[8]
	Ke = params[9]
	dG1 = params[10]
	kdeg = params[11]

	derivs = [
	- 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A)))),
	+ 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A)))) - 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))),
	+ 1.0 * (kcat1*E*(A/Ka)/(1+A/Ka+B/Kb*(M/Km))*(1-exp(dG0/(R*T)+log(B*M/A)))) - 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))),
	+ 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) - 1.0 * (kdeg*E),
	+ 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) - 1.0 * (kdeg*R),
	+ 1.0 * (kcat2*R*(B/Kb)*(M/Km)/(1+B/Kb*(M/Km)+R/Kr*(W/Kw)*(E/Ke))*(1-exp(dG1/(R*T)+log(R*W*E/(B*M))))) + 1.0 * (kdeg*R) + 1.0 * (kdeg*E)]
	return derivs

keysVar = ['A','B','M','E','R','W']
valuesVar = [50.0,0.1,0.1,0.05,8.3145,0.1]
dictVar = dict(zip(keysVar,valuesVar))

keysPar = ['kcat1','Ka','Kb','Km','dG0','T','kcat2','Kr','Kw','Ke','dG1','kdeg']
valuesPar = [10.0,0.1,0.1,0.1,-1000.0,298.0,10.0,0.1,0.1,0.1,-1000.0,0.1]
dictPar = dict(zip(keysPar,valuesPar))

