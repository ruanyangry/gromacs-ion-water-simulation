# _*_ coding:utf-8 _*_

import Gmxmdp as GMP
import FF 
import os

amber03={"IB+":1.000,"CA":2.000,"CL":-1.000,"NA":1.000,"MG":2.000,\
"K":1.000,"RB":1.000,"CS":1.000,"LI":1.000,"ZN":2.000,"water":["spce",\
"spc","tip3p","tip4p","tip4pew","tip5p"]}

charmm27={"NA":1.000,"MG":2.000,"K":1.000,"CS":1.000,"CA":2.000,"CL":-1.000,\
"ZN":2.000,"water":["spc","spce","tip3p","tip4p","tip5p","tips3p"]}

charmm36={"LI":[1.000,"LIT","LI"],"NA":[1.000,"SOD","NA"],"K":[1.000,"POT","K"],\
"CS":[1.000,"CES","Ces"],"CL":[-1.000,"CLA","CL"],"CA":[2.000,"CAL","Cal"],\
"MG":[2.000,"MGA","MG"],"ZN":[2.000,"ZN","ZN"],"water":["spc","spce","tip3p","tip4p"]}

oplsaa={"MG":2.000,"CA":2.000,"LI":1.000,"NA":1.000,"K":1.000,"RB":1.000,\
"CS":1.000,"F":-1.000,"CL":-1.000,"BR":-1.000,"I":-1.000,"water":["spc",\
"spce","tip3p","tip4p","tip4pew","tip5p","tip5pe"]}

gromos54a7={"CU1":1.000,"CU":2.000,"ZN":2.000,"MG":2.000,"CA":2.000,"NA":1.000,\
"CL":-1.000,"water":["spc","spce","tip3p","tip4p"]}

FF.writetop("amber03.ff",amber03["water"][2],["NA",amber03["NA"],1],\
["CL",amber03["CL"],1],16)

GMP.nemd(ensemble='nvt',integrator='md',nsteps=500000,dt=0.002,outfrequency=500,\
constraints='all-bonds',rcoulomb=1.2,rvdw=1.2,tcoupl='V-rescale',ref_t=300.0,pcoupl='Parrinello-Rahman',\
pcoupltype='isotropic',ref_p=1.0,nemd=["electric-field",1.0])

GMP.Equilibrium(ensemble='nvt',integrator='md',nsteps=500000,dt=0.002,outfrequency=500,\
constraints='all-bonds',rcoulomb=1.2,rvdw=1.2,tcoupl='V-rescale',ref_t=300.0,pcoupl='no',\
pcoupltype='isotropic',ref_p=1.0)

GMP.Equilibrium(ensemble='npt',integrator='md',nsteps=2000000,dt=0.002,outfrequency=500,\
constraints='all-bonds',rcoulomb=1.2,rvdw=1.2,tcoupl='V-rescale',ref_t=300.0,pcoupl='Parrinello-Rahman',\
pcoupltype='isotropic',ref_p=1.0)

os.system("dos2unix ionwater.sh")
os.system("qsub ionwater.sh")
