# _*_ coding:utf-8 _*_

'''
Author: Ruan Yang
Email: ruanyang_njut@163.com
Date: 2018.6.3
'''

def makendx(filename,cationname,anionname):
	'''
	Generated .ndx file for gromacs analysis tools.
	
	filename: the name of the .gro.
	cationname: the name of the cation.
	anionname: the name of the anion.
	'''
	with open(filename,"r") as f:
		lines=f.readlines()
		natoms=int(lines[1])
	with open("ndx.ndx","w") as f:
		# system
		f.write("[ system ]\n")
		system=[i+1 for i in range(natoms)]
		j=0
		for i in range(len(system)):
			j+=1
			print("%d"%(system[i]),end=' ',file=f)
			if j%20==0:
				print(end="\n",file=f)
		f.write("\n")
		# OW
		OW=[]
		HW=[]
		cation=[]
		anion=[]
		for line in lines:
			words=line.strip().split()
			if len(words) > 2:
				if words[1]=="OW":
					OW.append(int(words[2]))
				if words[1]=="HW1" or words[1]=="HW2":
					HW.append(int(words[2]))
				if words[1]==cationname:
					cation.append(int(words[2]))
				if words[1]==anionname:
					anion.append(int(words[2]))
				
		j=0
		f.write("[ OW ]\n")
		for i in range(len(OW)):
			j+=1
			print("%d"%(OW[i]),end=' ',file=f)
			if j%20==0:
				print(end="\n",file=f)
		f.write("\n")
		
		j=0
		f.write("[ HW ]\n")
		for i in range(len(HW)):
			j+=1
			print("%d"%(HW[i]),end=' ',file=f)
			if j%20==0:
				print(end="\n",file=f)
		f.write("\n")
		
		f.write("[ %s ]\n"%(cationname))
		for i in range(len(cation)):
			print("%d"%(cation[i]),end=' ',file=f)
		f.write("\n")
		
		f.write("[ %s ]\n"%(anionname))
		for i in range(len(anion)):
			print("%d"%(anion[i]),end=' ',file=f)
		f.write("\n")
		
		groups={}
		groups["system"]=0
		groups["OW"]=1
		groups["HW"]=2
		groups["%s"%(cationname)]=3
		groups["%s"%(anionname)]=4
		return groups,OW,HW,cation,anion
				
#groups,OW,HW,cation,anion=makendx("nvtnemd.gro","NA","CL")
#print(groups)
