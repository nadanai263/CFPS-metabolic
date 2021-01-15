

# Reaction 1
# 1.5 atp + 1.5 glc -> 1.5 g6p + 1.5 adp, catalyzed by Hex
# variables atp, glc, adp, g6p
# parameters VmaxHex, Katp1, Kglc1, Kadp1, Kg6p1, Keq1
VmaxHex * (atp*glc - g6p*adp/Keq1)/(Katp1*Kglc1) * 1/((1+atp/Katp1)*(1+glc/Kglc1)+(1+adp/Kadp1)*(1+g6p/Kg6p1)-1)

# Reaction 2
# 1.5 g6p -> 1.5 f6p, catalyzed by Pgi
# variables g6p, f6p
# parameters VmaxPgi, Kg6p2, Kf6p2, Keq2
VmaxPgi * (g6p - f6p/Keq2)/(Kg6p2) * 1/((1+g6p/Kg6p2)+(1+f6p/Kf6p2)-1)

# Reaction 3
# 1.5 f6p + 1.5 atp -> 1.5 fbp + 1.5 adp, catalyzed by Pfk
# variables f6p, atp, fbp, adp
# parameters VmaxPfk, Kadp3, Kf6p3, Kadp3, Kfbp3, Keq3
VmaxPfk * (f6p*atp - fbp*adp/Keq3)/(Katp3*Kf6p3) * 1/((1+atp/Katp3)*(1+f6p/Kf6p3)+(1+adp/Kadp3)*(1+fbp/Kfbp3)-1)

# Reaction 4
# 1.5 fbp -> 3 gap + 1.5 dhap, catalyzed by Fba
# variables fbp, gap, dhap
# parameters VmaxFba, Kfbp4, Kdhap4, Kgap4, Keq4
VmaxFba * (fbp - gap*dhap/Keq4)/(Kfbp4) * 1/((1+fbp/Kfbp4)+(1+dhap/Kdhap4)*(1+gap/Kgap4)-1)

# Reaction 5
# 1.5 dhap -> 3 gap, catalyzed by Tpi
# variables dhap, gap
# parameters VmaxTpi, Kdhap5, Kgap5, Keq5
VmaxTpi * (dhap - gap/Keq5)/Kdhap5 * 1/((1+dhap/Kdhap5)+(1+gap/Kgap5)-1)

# Reaction 6
# 3 gap + 3 ppi + nad -> 3 bpg + nadh, catalyzed by Gap
# variables gap, nad, ppi, bpg, nadh
# parameters VmaxGap, Kgap6, Knad6, Kppi6, Kbpg6, Knadh6, Keq6
VmaxGap * (gap*nad*ppi - bpg*nadh/Keq6)/(Kgap6*Knad6*Kppi6) * 1/((1+gap/Kgap6)*(1+nad/Knad6)*(1+ppi/Kppi6)+bpg/Kbpg6*(1+nadh/Knadh6)-1)

# Reaction 7
# 3 gap + 3 ppi + 2 nadp -> 3 bpg + 2 nadph, catalyzed by mGap
# variables gap, nadp, ppi, bpg, nadph
# parameters VmaxGap, Kgap7, Knadp7, Kppi7, Kbpg7, Knadph7, Keq7
VmaxmGap * (gap*nadp*ppi - bpg*nadph/Keq7)/(Kgap7*Knadp7*Kppi7) * 1/((1+gap/Kgap7)*(1+nadp/Knadp7)*(1+ppi/Kppi7)+bpg/Kbpg7*(1+nadph/Knadph7)-1)

# Reaction 8
# 3 bpg + 3 adp -> 3 p3g + 3 atp, catalyzed by Pgk
# variables bpg, adp, p3g, atp
# parameters VmaxPgk, Kbpg8, Kadp8, Kp3g8, Katp8, Keq8
VmaxPgk * (bpg*adp - p3g*atp/Keq8)/(Kbpg8*Kadp8) * 1/((1+bpg/Kbpg8)*(1+adp/Kadp8)+(1+p3g/Kp3g8)*(1+atp/Katp8)-1)

# Reaction 9
# 3 p3g -> 3 p2g, catalyzed by Pgm
# variables p3g, p2g
# parameters VmaxPgm, Kp3g9, Kp2g9, Keq9
VmaxPgm * (p3g - p2g/Keq9)/Kp3g9 * 1/((1+p3g/Kp3g9)+(1+p2g/Kp2g9)-1)

# Reaction 10
# 3 p2g -> 3 pep, catalyzed by Eno
# variables p2g, pep
# parameters VmaxEno, Kp2g10, Kpep10, Keq10
VmaxEno * (p2g - pep/Keq10)/Kp2g10 * 1/((1+p2g/Kp2g10)+(1+pep/Kpep10)-1)

# Reaction 11
# 3 pep + 3 adp -> 3 pyr + 3 atp, catalyzed by Pyk
# variables pep, adp, pyr, atp
# parameters VmaxPyk, Kadp11, Kpep11, Katp11, Kpyr11, Keq11
VmaxPyk * (pep*adp - pyr*atp/Keq11)/(Kadp11*Kpep11) * 1/((1+adp/Kadp11)*(1+pep/Kpep11)+(1+atp/Katp11)*(1+pyr/Kpyr11)-1)

# Reaction 12
# 3 pyr + 3 coa + 3 nad -> 3 accoa + 3 nadh + 3 co2, catalyzed by Pdh
# variables coa, nad, pyr, accoa, co2, nadh 
# parameters VmaxPdh, Kcoa12, Knad12, Kpyr12, Kaccoa12, Kco212, Knadh12, Keq12
VmaxPdh * (coa*nad*pyr - accoa*co2*nadh/Keq12)/(Kcoa12*Knad12*Kpyr12) * 1/((1+coa/Kcoa12)*(1+nad/Knad12)*(1+pyr/Kpyr12)+(1+accoa/Kaccoa12)*(1+co2/Kco212)*(1+nadh/Knadh12)-1)

# Reaction 13 ERROR IN SI FORMULA?
# 2 accoa -> coa + acaccoa, catalyzed by Thl (PhaA)
# variables accoa, acaccoa, coa
# parameters VmaxThl, Kaccoa13, Kacaccoa13, Kcoa13, Keq13
VmaxThl * (accoa*accoa - coa*acaccoa/Keq13)/(Kaccoa13*Kaccoa13) * 1/((1+accoa/Kaccoa13)*(1+accoa/Kaccoa13)+(1+acaccoa/Kacaccoa13)*(1+coa/Kcoa13)-1)

# Reaction 14
# acaccoa + accoa -> coa + hmgcoa, catalyzed by Hmgs
# variables acaccoa, accoa, coa, hmgcoa
# parameters VmaxHmgs, Kacaccoa14, Kaccoa14, Khmgcoa14, Kcoa14, Keq14
VmaxHmgs * (acaccoa*accoa - coa*hmgcoa/Keq14)/(Kacaccoa14*Kaccoa14) * 1/((1+acaccoa/Kacaccoa14)*(1+accoa/Kaccoa14)+(1+hmgcoa/Khmgcoa14)*(1+coa/Kcoa14)-1)

# Reaction 15
# hmgcoa + 2 nadph -> mev + coa + 2 nadp, catalyzed by Hmgr
# variables hmgcoa, nadph, mev, coa, nadp
# parameters VmaxHmgr, Khmgcoa15, Knadph15, Kmev15, Kcoa15, Knadp15, Keq15
VmaxHmgr * (hmgcoa*nadph*nadph - mev*coa*nadp*nadp/Keq15)/(Khmgcoa15*Knadph15*Knadph15) * 1/((1+hmgcoa/Khmgcoa15)*(1+nadph/Knadph15)*(1+nadph/Knadph15)+(1+mev/Kmev15)*(1+coa/Kcoa15)*(1+nadp/Knadp15)*(1+nadp/Knadp15)-1)

# Reaction 16
# mev + atp -> mvp + adp, catalyzed by Mvk
# variables mev, atp, mvp, adp
# parameters VmaxMvk, Kmev16, Katp16, Kmvp16, Kadp16, Keq16
VmaxMvk * (mev*atp - mvp*adp/Keq16)/(Kmev16*Katp16) * 1/((1+mev/Kmev16+mvp/Kmvp16)*(1+atp/Katp16+adp/Kadp16))

# Reaction 17
# mvp + atp -> mpp + adp, catalyzed by Pmvk
# variables mvp, atp, mpp, adp
# parameters VmaxPmvk, Kmvp17, Katp17, Kmpp17, Kadp17, Keq17
VmaxPmvk * (mvp*atp - mpp*adp/Keq17)/(Kmvp17*Katp17) * 1/((1+mvp/Kmvp17+mpp/Kmpp17)*(1+atp/Katp17+adp/Kadp17))

# Reaction 18
# mpp + atp -> ipp + co2 + adp + ppi, catalyzed by Mdc
# variables mpp, atp, adp, co2, ipp, ppi
# parameters VmaxMdc, Kmpp18, Katp18, Kadp18, Kco218, Kipp18, Kppi18, Keq18
VmaxMdc * (mpp*atp - adp*co2*ipp*ppi/Keq18)/(Kmpp18*Katp18) * 1/((1+mpp/Kmpp18)*(1+atp/Katp18)+(1+adp/Kadp18)*(1+co2/Kco218)*(1+ipp/Kipp18)*(1+ppi/Kppi18)-1)

# Reaction 19
# ipp -> dmapp catalyzed by Idi
# variables ipp, dmapp
# parameters VmaxIdi, Kipp19, Kdmapp19, Keq19
VmaxIdi * (ipp - dmapp/Keq19)/Kipp19 * 1/((1+ipp/Kipp19)+(1+dmapp/Kdmapp19)-1)

# Reaction 20
# dmapp + ipp -> gpp + pppi, catalyzed by Fpps-S82F
# variables dmapp, ipp, gpp, pppi
# parameters VmaxFS82F, Kdmapp20, Kipp20, Kgpp20, Kpppi20, Keq20
VmaxFS82F * (dmapp*ipp - gpp*pppi/Keq20)/(Kdmapp20*Kipp20) * 1/((1+dmapp/Kdmapp20)*(1+ipp/Kipp20)+(1+gpp/Kgpp20)*(1+pppi/Kpppi20)-1)

# Reaction 21
# gpp -> lim + pppi, catalyzed by Lim
# variables gpp, lim, pppi
# parameters VmaxLim, Kgpp21, Klim21, Kpppi21, Keq21
VmaxLim * (gpp - lim*pppi/Keq21)/Kgpp21 * 1/((1+gpp/Kgpp21)+(1+lim/Klim21)*(1+pppi/Kpppi21)-1)

# Reaction 22
# pppi -> 2 ppi, catalyzed by Ppase
# variables pppi, ppi
# parameters VmaxPpase, Kpppi22, Kppi22, Keq22
VmaxPpase * (pppi - ppi*ppi/Keq22)/Kpppi * 1/((1+pppi/Kpppi22)+(1+ppi/Kppi22)*(1+ppi/Kppi22)-1)

# Reaction 23
# 4 nadh -> 4 nad, catalyzed by NoxE
# variables nadh, nad
# parameters VmaxNoxe, Knadh23
VmaxNoxe * nadh/(Knadh23+nadh)