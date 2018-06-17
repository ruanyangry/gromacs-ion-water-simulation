# gromacs-ion-water-simulation
This repository contained  python code used to do ion water system simulation.

# Authors:  

Ruan Yang  
Email: ruanyang_njut@163.com  

Required library:  

MDAnalysis； https://github.com/MDAnalysis/mdanalysis  

# Uage:  

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

# Analysis methods in ionwateranalysis-v3.py:  

1. Radial distribution functions for ow-ow、ow-hw、cation-ow and anion-ow. 
2. The coordination numbers for the first and second shells of ions.  
3. Probability distributions of water dipole orientation.-OH and -HH.  
4. Probability distributions of the orientation angle α and θ of water molecules in the first and second hydration shells.  
5. Analyzing hydrogen bond lifetimes (HBL) and hydrogen bond (HB) in ions first or second hydration shells and between first and second hydration shells.  
6. Mean square displacement (MSD) of ions and water along x y z and the total displacment.  
7. Analyzing survival probability (SP) SurvivalProbability for water molecules in ions hydration shells.  

The analysis codes can automatically obtain the radius of the ions first and second hydration layer.  

# Citation:  
These publications associated with this code is found here:

1. Molecular Dynamics Study of Mg2+/Li+ Separation via Biomimetic Graphene-Based Nanopores: The Role of Dehydration in Second Shell. DOI: 10.1021/acs.langmuir.6b03001.  

2. Mg2+-Channel-Inspired Nanopores for Mg2+/Li+ Separation: The Effect of Coordination on the Ionic Hydration Microstructure. DOI: 10.1021/acs.langmuir.7b01249.  



