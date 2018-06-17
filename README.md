# gromacs-ion-water-simulation
This repository contained  python code used to do ion water system simulation.

# Authors:  

Ruan Yang  
Email: ruanyang_njut@163.com  

Required library:  

MDAnalysis； https://github.com/MDAnalysis/mdanalysis  

# Uage:  
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
    
    Gmxmdp.py: Generated GROMACS .mdp file. Simulation process: EM---NVTEquilibrium---NPTEquilibrium---NVT(NPT) Production run---Standmdp
    
    FF.py: Generated GROMACS .top file and ionwater.sh.

# Analysis methods in ionwateranalysis-v3.py:  

1. Radial distribution functions for ow-ow、ow-hw、cation-ow and anion-ow. 
2. The coordination numbers for the first and second shells of ions.  
3. Probability distributions of water dipole orientation.-OH and -HH.  
4. Probability distributions of the orientation angle α and θ of water molecules in the first and second hydration shells.  
5. Analyzing hydrogen bond lifetimes (HBL) and hydrogen bond (HB) in ions first or second hydration shells and between first and second hydration shells.  
6. Mean square displacement (MSD) of ions and water along x y z and the total displacment.  
7. Analyzing survival probability (SP) SurvivalProbability for water molecules in ions hydration shells.  

The analysis codes can automatically obtain the radius of the ions first and second hydration layer.  

    python ry.py -gro "nvtnemd.gro" -trr "nvtnemd.trr" -tpr "nvtnemd.tpr" -cationname "NA" -anionname "CL"  

    parser=argparse.ArgumentParser(description="Analysis ion water systems.\
    Using MDAnalysis library.")
    parser.add_argument('-gro',dest='GRO',type=str,help='The name of .gro file')
    parser.add_argument('-trr',dest='TRR',type=str,help='The name of .trr file')
    parser.add_argument('-tpr',dest='TPR',type=str,help='The name of .tpr file')
    parser.add_argument('-cationname',dest='cationname',type=str,help='the name of cation ions')
    parser.add_argument('-anionname',dest='anionname',type=str,help='the name of anion ions')
    
    args = parser.parse_args()
    GRO=args.GRO
    TRR=args.TRR
    TPR=args.TPR
    cationname=args.cationname
    anionname=args.anionname
    
References: https://www.mdanalysis.org/docs/documentation_pages/analysis/waterdynamics.html 

# Citation:  
These publications associated with this code is found here:

1. Molecular Dynamics Study of Mg2+/Li+ Separation via Biomimetic Graphene-Based Nanopores: The Role of Dehydration in Second Shell. DOI: 10.1021/acs.langmuir.6b03001.  

2. Mg2+-Channel-Inspired Nanopores for Mg2+/Li+ Separation: The Effect of Coordination on the Ionic Hydration Microstructure. DOI: 10.1021/acs.langmuir.7b01249.  



