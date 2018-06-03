# _*_ coding:utf-8 _*_

'''
Author: Ruan Yang
Email: ruanyang_njut@163.com
'''

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

def writetop(ff,watertype,cation,anion,nt):
	'''
	Use gromacs commands gmx solvate and gmx genion	
	
	ff: "amber03.ff","charmm27.ff","ccharmm36.ff","oplsaa.ff","gromos54a7.ff"
	watertype: select water types
	cation(list): select cation types
	              cation[0]:atom name
	              cation[1]:atom charge
	              cation[2]:atom number
	anion(list): select anion types
	              anion[0]:atom name
	              anion[1]:atom charge
	              anion[2]:atom number
	nt: default=16
	'''
	with open("water.top","w") as f:
		f.write('; Author: Ruan Yang\n')
		f.write('; Email: ruanyang_njut@163.com\n')
		f.write('; This top file for ion-water system.\n')
		f.write('\n')
		f.write('#include "%s/forcefield.itp"\n'%(ff))
		f.write('\n')
		f.write('#include "%s/%s.itp"\n'%(ff,watertype))
		f.write("\n")
		f.write('#include "%s/%s.itp"\n'%(ff,"ions"))
		f.write("\n")
		f.write('[ system ]\n')
		f.write('   Water-Ion\n')
		f.write('\n')
		f.write('[ molecules ]\n')
		f.write('; Compound        nmols\n')
		f.write('   SOL       2138\n')
	with open("emptybox.gro","w") as f:
		f.write("Empty box\n")
		f.write("  0\n")
		f.write("  4.0  4.0 4.0\n")
	with open("em.mdp","w") as f:
		f.write('; Author: Ruan Yang\n')
		f.write('; Email: ruanyang_njut@163.com\n')
		f.write("; Used as input into grompp to generate em.tpr\n")
		f.write("\n")
		f.write("integrator  =  steep\n")
		f.write("emtol       =  1000.0\n")
		f.write("emstep      =  0.01\n")
		f.write("nsteps      =  50000\n")
		f.write("\n")
		f.write("nstlist     =  1\n")
		f.write("cutoff-scheme = Verlet\n")
		f.write("ns_type     =  grid\n")
		f.write("coulombtype =  PME\n")
		f.write("rcoulomb    =  1.0\n")
		f.write("rvdw        =  1.0\n")
		f.write("pbc         =  xyz\n")
		
	#import os
	## Use gromacs commands gmx solvate and gmx genion	
	#os.system("gmx solvate -cp emptybox.gro -cs -o water.gro")
	#os.system("gmx grompp -f em.mdp -c water.gro -p water.top -o em.tpr")
	#os.system("echo 2 | gmx genion -s em.tpr -o waterion.gro -pname %s \
	#-pq %d -np %d -nname %s -nq %d -nn %d -p water.top"%(cation[0],cation[1],\
	#cation[2],anion[0],anion[1],anion[2]))
	#os.system("gmx grompp -f em.mdp -c waterion.gro -p water.top -o em.tpr")
	#os.system("gmx mdrun -v -deffnm em")
	#os.system("gmx grompp -f nvtequilibrium.mdp -c em.gro -p water.top -o \
	#nvtequilibrium.tpr -maxwarn 2")
	#os.system("gmx mdrun -v -deffnm nvtequilibrium")
	#os.system("gmx grompp -f nptequilibrium.mdp -c nvtequilibrium.gro -p water.top -o \
	#nptequilibrium.tpr -maxwarn 2")
	#os.system("gmx mdrun -v -deffnm nptequilibrium")
	#os.system("gmx grompp -f nvtnemd.mdp -c nptequilibrium.gro -p water.top -o \
	#nvtnemd.tpr -maxwarn 2")
	#os.system("gmx mdrun -v -deffnm nvtnemd")
	with open("ionwater.sh","w") as f:
		f.write("#!/bin/bash\n")
		f.write('#$ -S /bin/sh\n')
		f.write('#$ -N test\n')
		f.write('#$ -j y\n')
		f.write('#$ -o .\n')
		f.write('#$ -e .\n')
		f.write('#$ -cwd\n')
		f.write('#$ -q all.q\n')
		f.write('#$ -masterq all.q@node9\n')
		f.write('#$ -pe thread %d-%d\n'%(nt,nt))
		f.write('source ~/.bashrc\n')
		f.write('hash -r\n')
		#f.write('export PATH=/export/home/ry/gromacs-5.0/bin:PATH\n')
		f.write('\n')
		f.write("\n")
		f.write("# Author: Ruan Yang\n")
		f.write("# Email: ruanyang_njut@163.com\n")
		f.write("\n")
		f.write("mkdir model\n")
		f.write("cd model\n")
		f.write("gmx solvate -cp ../emptybox.gro -cs -o water.gro\n")
		f.write("wait\n")
		f.write("\n")
		f.write("gmx grompp -f ../em.mdp -c water.gro -p ../water.top -o em.tpr\n")
		f.write("wait\n")
		f.write("\n")
		f.write("echo 2 | gmx genion -s em.tpr -o waterion.gro -pname %s \
		-pq %d -np %d -nname %s -nq %d -nn %d -p ../water.top\n"%(cation[0],\
		cation[1],cation[2],anion[0],anion[1],anion[2]))
		f.write("wait\n")
		f.write("cd ../\n")
		f.write("\n")
		f.write("mkdir em\n")
		f.write("cd em\n")
		f.write("gmx grompp -f ../em.mdp -c ../model/waterion.gro -p ../water.top -o em.tpr\n")
		f.write("gmx mdrun -v -nt $NSLOTS -deffnm em\n")
		f.write("wait\n")
		f.write("cd ../\n")
		f.write("\n")
		f.write("mkdir nvte\n")
		f.write("cd nvte\n")
		f.write("gmx grompp -f ../nvtequilibrium.mdp -c ../em/em.gro -p ../water.top -o \
		nvtequilibrium.tpr -maxwarn 2\n")
		f.write("gmx mdrun -v -nt $NSLOTS -deffnm nvtequilibrium\n")
		f.write("wait\n")
		f.write("cd ../\n")
		f.write("\n")
		f.write("mkdir npte\n")
		f.write("cd npte\n")
		f.write("gmx grompp -f ../nptequilibrium.mdp -c ../nvte/nvtequilibrium.gro -p \
		../water.top -t ../nvte/nvtequilibrium.cpt -o nptequilibrium.tpr -maxwarn 2\n")
		f.write("gmx mdrun -v -nt $NSLOTS -deffnm nptequilibrium\n")
		f.write("wait\n")
		f.write("cd ../\n")
		f.write("\n")
		f.write("mkdir nemd\n")
		f.write("cd nemd\n")
		f.write("gmx grompp -f ../nvtnemd.mdp -c ../npte/nptequilibrium.gro -p ../water.top \
		-t ../npte/nptequilibrium.cpt -o nvtnemd.tpr -maxwarn 2\n")
		f.write("gmx mdrun -v -nt $NSLOTS -deffnm nvtnemd\n")
		f.write("cd ../\n")
		f.write("wait\n")		

#writetop("amber03.ff",amber03["water"][2],["NA",amber03["NA"],1],\
#["CL",amber03["CL"],1],16)
#
#import os
#os.system("dos2unix ionwater.sh")
#os.system("qsub ionwater.sh")
