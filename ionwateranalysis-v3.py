# _*_ coding:utf-8 _*_

'''
Author: Ruan Yang
Email: ruanyang_njut@163.com
References: https://www.mdanalysis.org/
            https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html
            https://www.mdanalysis.org/docs/documentation_pages/analysis/rdf.html
            https://www.mdanalysis.org/docs/documentation_pages/analysis/waterdynamics.html
            https://www.mdanalysis.org/docs/documentation_pages/analysis/lineardensity.html
Just for GROMACS topology and trajectory
'''

import MDAnalysis as mda
mda.core.flags['use_periodic_selections'] = True
from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR
from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP
from MDAnalysis.analysis.rdf import InterRDF as IRDF
from MDAnalysis.analysis.lineardensity import LinearDensity as LA
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import collections
import os

print("#------------------------------------------------------#")
print("Import module done")
print("#------------------------------------------------------#")
print("\n")

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

# Defined variable value

cation=[]
cation.append(cationname)
anion=[]
anion.append(anionname)
u=mda.Universe(GRO,TRR)
start=0
end=u.trajectory.n_frames

# RDF analysis

print("#------------------------------------------------------#")
print("Start RDF analysis")
os.system("mkdir rdf")
os.system("cd rdf")

ow=u.select_atoms("name OW")
hw=u.select_atoms("name HW1 or name HW2")
selcation=u.select_atoms("name %s"%(cation[0]))
selanion=u.select_atoms("name %s"%(anion[0]))

ow_ow=IRDF(ow,ow,nbins=100,exclusion_block=(1,1))
rdf_ow_ow=ow_ow.run()
print("ow ow rdf analysis done")
ow_hw=IRDF(ow,hw,nbins=100)
rdf_ow_hw=ow_hw.run()
print("ow hw rdf analysis done")
cation_ow=IRDF(selcation,ow,nbins=100)
rdf_cation_ow=cation_ow.run()
print("cation ow rdf analysis done")
anion_ow=IRDF(selanion,ow,nbins=100)
rdf_anion_ow=anion_ow.run()
print("anion ow rdf analysis done")

with open("rdf.txt","w") as f:
	f.write("# RDF result for ow-ow,ow-hw,cation-ow,anion-ow\n")
	for i in range(len(ow_ow.bins)):
		f.write("%.3f %.3f %.3f %.3f %.3f \n"%(rdf_ow_ow.bins[i],\
		rdf_ow_ow.rdf[i],rdf_ow_hw.rdf[i],rdf_cation_ow.rdf[i],rdf_anion_ow.rdf[i]))
		
data=np.loadtxt("rdf.txt",skiprows=1)	

max_index_f=np.argmax(data[:,-2])
min_index_f=np.argmin(data[max_index_f:,-2])
cationfirst=data[max_index_f+min_index_f,0]
cation.append(cationfirst)

max_index_s=np.argmax(data[max_index_f+min_index_f+1:,-2])
min_index_s=np.argmax(data[max_index_f+min_index_f+1+max_index_s:,-2])
cationsecond=data[max_index_f+min_index_f+1+max_index_s+min_index_s,0]
cation.append(cationsecond)

print("%s first hydration radius = %.4f"%(cation[0],cation[1]))
print("%s second hydration radius = %.4f"%(cation[0],cation[2]))

max_index_f=np.argmax(data[:,-1])
min_index_f=np.argmin(data[max_index_f:,-1])
anionfirst=data[max_index_f+min_index_f,0]
anion.append(anionfirst)

max_index_s=np.argmax(data[max_index_f+min_index_f+1:,-1])
min_index_s=np.argmax(data[max_index_f+min_index_f+1+max_index_s:,-1])
anionsecond=data[max_index_f+min_index_f+1+max_index_s+min_index_s,0]
anion.append(anionsecond)

print("%s first hydration radius = %.4f"%(anion[0],anion[1]))
print("%s second hydration radius = %.4f"%(anion[0],anion[2]))

print(cation)
print(anion)

plt.figure(1,figsize=(18,8))
plt.subplot(141)
plt.plot(data[:,0],data[:,1],'g-',lw=2.0,label="ow-ow rdf")
plt.xlabel("distance")
plt.ylabel("RDF")
plt.title("ow-ow rdf")

plt.subplot(142)
plt.plot(data[:,0],data[:,2],'g-',lw=2.0,label="ow-hw rdf")
plt.xlabel("distance")
plt.ylabel("RDF")
plt.title("ow-hw rdf")

plt.subplot(143)
plt.plot(data[:,0],data[:,3],'g-',lw=2.0,label="%s-ow rdf"%(cation[0]))
plt.xlabel("distance")
plt.ylabel("RDF")
plt.title("%s-ow rdf"%(cation[0]))

plt.subplot(144)
plt.plot(data[:,0],data[:,4],'g-',lw=2.0,label="%s-ow rdf"%(anion[0]))
plt.xlabel("distance")
plt.ylabel("RDF")
plt.title("%s-ow rdf"%(anion[0]))

plt.savefig("rdf.jpg",dpi=300)
plt.clf()

os.system("mv *.txt *.jpg rdf")

print("RDF analysis done")
print("#------------------------------------------------------#")
print("\n")

sel=[]
sel.append("name OW and (around %.4f resname %s)"%(cation[1],cation[0]))
sel.append("name OW and (around %.4f resname %s) and name OW and (not \
around %.4f resname %s)"%(cation[2],cation[0],cation[1],cation[0]))
sel.append("name OW and (around %.4f resname %s)"%(anion[1],anion[0]))
sel.append("name OW and (around %.4f resname %s) and name OW and (not \
around %.4f resname %s)"%(anion[2],anion[0],anion[1],anion[0]))

# Water mass density analysis

print("#------------------------------------------------------#")
print("Start Water mass density analysis")
os.system("mkdir massdensity")

resid=int((u.select_atoms("resname SOL").n_atoms)/3)
u=mda.Universe(TPR,TRR)
water=u.select_atoms("resid 0-%d"%(resid))
watermassdensity=LA(water)
watermassdensity.run()
watermassdensity.save(description='watermassdensity', form='txt')

data=np.loadtxt("nvtnemd.watermassdensity_atoms.ldens",skiprows=2)

labels=["x","y","z"]
colors=["chartreuse","deepskyblue","r"]
for i in range(len(labels)):
	if i == 0:
		plt.plot(data[:,0]/10.0,data[:,i+1],color=colors[i],lw=2.0,label=labels[i])
	if i == 1:
		plt.plot(data[:,0]/10.0,data[:,i+2],color=colors[i],lw=2.0,label=labels[i])
	if i ==2:
		plt.plot(data[:,0]/10.0,data[:,i+3],color=colors[i],lw=2.0,label=labels[i])
plt.xlabel("Distance (nm)")
plt.ylabel("Mass Density (g/cm^3)")
plt.legend(loc="best")
plt.title("Water mass density along x y z axis")
plt.savefig("Water-mass-density.jpg",dpi=300)
plt.clf()

for i in range(len(labels)):
	if i == 0:
		plt.plot(data[:,0]/10.0,data[:,i+7],color=colors[i],lw=2.0,label=labels[i])
	if i == 1:
		plt.plot(data[:,0]/10.0,data[:,i+8],color=colors[i],lw=2.0,label=labels[i])
	if i == 2:
		plt.plot(data[:,0]/10.0,data[:,i+9],color=colors[i],lw=2.0,label=labels[i])
plt.xlabel("Distance (nm)")
plt.ylabel("Charge Density (e/A^3)")
plt.legend(loc="best")
plt.title("Water charge density along x y z axis")
plt.savefig("Water-charge-density.jpg",dpi=300)

os.system("mv *.ldens *.jpg massdensity")
plt.clf()

print("Water mass density analysis done")
print("#------------------------------------------------------#")
print("\n")

u=mda.Universe(GRO,TRR)

# Hydration number analysis

print("#------------------------------------------------------#")
print("Start Hydration number analysis")
os.system("mkdir Hydration-number")

hn=np.zeros([end,5])
cationnumber=u.select_atoms("resname %s"%(cation[0])).n_atoms
anionnumber=u.select_atoms("resname %s"%(anion[0])).n_atoms

for ts in u.trajectory:
	hn[ts.frame,0]=ts.frame
	for i in range(len(sel)):
		hn[ts.frame,i+1]=u.select_atoms("%s"%(sel[i]),updating=True).n_atoms
	hn[ts.frame,1]=hn[ts.frame,1]/cationnumber
	hn[ts.frame,2]=hn[ts.frame,2]/cationnumber
	hn[ts.frame,3]=hn[ts.frame,3]/anionnumber
	hn[ts.frame,4]=hn[ts.frame,4]/anionnumber

np.savetxt("hydration-number.txt",hn,fmt="%d %.4f %.4f %.4f %.4f")

plt.figure(1,figsize=(18,6))
plt.subplot(121)
plt.plot(hn[:,0]/1000.0,hn[:,1],'g-',lw=2.0,label="%s first"%(cation[0]))
plt.plot(hn[:,0]/1000.0,hn[:,3],'b-',lw=2.0,label="%s first"%(anion[0]))
plt.xlabel("Time (ns)")
plt.ylabel("Hydration Number")
plt.ylim(min(np.min(hn[:,1]),np.min(hn[:,3]))-2.0,max(np.max(hn[:,1]),np.max(hn[:,3]))+2.0)
plt.legend(loc="best")
plt.title("Ions First hydration Number")

plt.subplot(122)
plt.plot(hn[:,0],hn[:,2],'g-',lw=2.0,label="%s second"%(cation[0]))
plt.plot(hn[:,0],hn[:,4],'b-',lw=2.0,label="%s second"%(anion[0]))
plt.xlabel("Time")
plt.ylabel("Hydration Number")
plt.ylim(min(np.min(hn[:,2]),np.min(hn[:,4]))-2.0,max(np.max(hn[:,2]),np.max(hn[:,4]))+2.0)
plt.legend(loc="best")
plt.title("Ions Second hydration Number")

plt.savefig("hydration-number.jpg",dpi=300)
plt.clf()

try:
	plt.figure(1,figsize=(18,8))
	pdf=collections.Counter(hn[:,1])
	results=[]
	for k,v in pdf.items():
		results.append([k,v])
	results=np.array(results)
	plt.bar(results[:,0],results[:,1],width=0.35)
	plt.xlabel("Hydration number")
	plt.ylabel("Occurrence number")
	averagehn=0.0
	for i in range(len(results)-1):
		averagehn += results[i+1,0]*results[i+1,1]
	averagehn=averagehn/sum(results[1:,1])
	plt.title("The average hydration number = %.4f"%(averagehn))
	plt.savefig("HN-PDF-%s-First.jpg"%(cation[0]),dpi=300)
	plt.clf()
	
	plt.figure(1,figsize=(18,8))
	pdf=collections.Counter(hn[:,3])
	results=[]
	for k,v in pdf.items():
		results.append([k,v])
	results=np.array(results)
	plt.bar(results[:,0],results[:,1],width=0.35)
	plt.xlabel("Hydration number")
	plt.ylabel("Occurrence number")
	averagehn=0.0
	for i in range(len(results)-1):
		averagehn += results[i+1,0]*results[i+1,1]
	averagehn=averagehn/sum(results[1:,1])
	plt.title("The average hydration number = %.4f"%(averagehn))
	plt.savefig("HN-PDF-%s-First.jpg"%(anion[0]),dpi=300)
	plt.clf()	
	
	plt.figure(1,figsize=(18,8))
	pdf=collections.Counter(hn[:,2])
	results=[]
	for k,v in pdf.items():
		results.append([k,v])
	results=np.array(results)
	plt.bar(results[:,0],results[:,1],width=0.35)
	plt.xlabel("Hydration number")
	plt.ylabel("Occurrence number")
	averagehn=0.0
	for i in range(len(results)-1):
		averagehn += results[i+1,0]*results[i+1,1]
	averagehn=averagehn/sum(results[1:,1])
	plt.title("The average hydration number = %.4f"%(averagehn))
	plt.savefig("HN-PDF-%s-Second.jpg"%(cation[0]),dpi=300)
	plt.clf()	
	
	plt.figure(1,figsize=(18,8))
	pdf=collections.Counter(hn[:,4])
	results=[]
	for k,v in pdf.items():
		results.append([k,v])
	results=np.array(results)
	plt.bar(results[:,0],results[:,1],width=0.35)
	plt.xlabel("Hydration number")
	plt.ylabel("Occurrence number")
	averagehn=0.0
	for i in range(len(results)-1):
		averagehn += results[i+1,0]*results[i+1,1]
	averagehn=averagehn/sum(results[1:,1])
	plt.title("The average hydration number = %.4f"%(averagehn))
	plt.savefig("HN-PDF-%s-Second.jpg"%(anion[0]),dpi=300)
	plt.clf()
except:
	print("Hydration number pdf analysis plot error")
	
os.system("mv *.txt *.jpg Hydration-number")

print("Hydration number analysis done")
print("#------------------------------------------------------#")
print("\n")

# hydrogen bond analysis

print("#------------------------------------------------------#")
print("Start hydrogen bond analysis")
os.system("mkdir hydrogen-bond")
os.system("cd hydrogen-bond")

sel=[]
sel.append("same resid as (resname SOL and (around %.4f resname %s))"%(cation[1],cation[0]))
sel.append("same resid as (resname SOL and (around %.4f resname %s) and \
resname SOL and (not around %.4f resname %s))"%(cation[2],cation[0],cation[1],cation[0]))
sel.append("same resid as (resname SOL and (around %.4f resname %s))"%(anion[1],anion[0]))
sel.append("same resid as (resname SOL and (around %.4f resname %s) and \
resname SOL and (not around %.4f resname %s))"%(anion[2],anion[0],anion[1],anion[0]))

for i in range(len(sel)):
	HB_analysis=HBA(u,sel[i],sel[i],distance=3.0, angle=120.0)
	HB_analysis.run()
	number=np.zeros([len(HB_analysis.count_by_time()),2])
	number[:,0]=np.array([column[0] for column in HB_analysis.count_by_time()]).\
	reshape(len(HB_analysis.count_by_time()))
	number[:,1]=np.array([column[1] for column in HB_analysis.count_by_time()]).\
	reshape(len(HB_analysis.count_by_time()))
	print("%s hydrogen bond analysis done"%(sel[i]))
	
	# Save Hbond:frame
	np.savetxt("hbond-frame-%d.txt"%(i),number,fmt="%10d %10d")
	
	try:
		plt.figure(1,figsize=(18,8))
		plt.plot(number[:,0],number[:,1],'g-',lw=2.0,label="HB")
		plt.xlabel("Time")
		plt.ylabel("HB Number")
		plt.title("Hydrogen bond number distribution")
		plt.legend(loc='best')
		plt.savefig("HB-Time-%d.jpg"%(i),dpi=300)
		plt.clf()
		
		plt.figure(1,figsize=(18,8))
		pdf=collections.Counter(number[:,1])
		results=[]
		for k,v in pdf.items():
			results.append([k,v])
		results=np.array(results)
		plt.bar(results[:,0],results[:,1],width=0.35)
		plt.xlabel("Hydrogen bond number")
		plt.ylabel("Occurrence number")
		averagehb=0.0
		for i in range(len(results)-1):
			averagehb += results[i+1,0]*results[i+1,1]
		averagehb=averagehb/sum(results[1:,1])
		plt.title("The average hydrogen bond = %.4f"%(averagehb))
		plt.savefig("HB-PDF-%d-%d.jpg"%(i,i+1),dpi=300)
		plt.clf()		
	except:
		print("Hydrogen bond in first or second hydration shell plot error")
	
for i in range(0,len(sel),2):
	HB_analysis=HBA(u,sel[i],sel[i+1],distance=3.0, angle=120.0)
	HB_analysis.run()
	number=np.zeros([len(HB_analysis.count_by_time()),2])
	number[:,0]=np.array([column[0] for column in HB_analysis.count_by_time()]).\
	reshape(len(HB_analysis.count_by_time()))
	number[:,1]=np.array([column[1] for column in HB_analysis.count_by_time()]).\
	reshape(len(HB_analysis.count_by_time()))
	print("%s %s hydrogen bond analysis done"%(sel[i],sel[i+1]))
	
	# Save Hbond:frame
	np.savetxt("hbond-frame-%d-%d.txt"%(i,i+1),number,fmt="%10d %10d")
	
	try:
		plt.figure(1,figsize=(18,8))
		plt.plot(number[:,0],number[:,1],'g-',lw=2.0,label="HB")
		plt.xlabel("Time")
		plt.ylabel("HB Number")
		plt.title("Hydrogen bond number distribution")
		plt.legend(loc='best')
		plt.savefig("HB-Time-%d-%d.jpg"%(i,i+1),dpi=300)
		plt.clf()
		
		
		plt.figure(1,figsize=(18,8))
		pdf=collections.Counter(number[:,1])
		results=[]
		for k,v in pdf.items():
			results.append([k,v])
		results=np.array(results)
		plt.bar(results[:,0],results[:,1],width=0.35)
		plt.xlabel("Hydrogen bond number")
		plt.ylabel("Occurrence number")
		averagehb=0.0
		for i in range(len(results)-1):
			averagehb += results[i+1,0]*results[i+1,1]
		averagehb=averagehb/sum(results[1:,1])
		plt.title("The average hydrogen bond = %.4f"%(averagehb))
		plt.savefig("HB-PDF-%d-%d.jpg"%(i,i+1),dpi=300)
		plt.clf()		
		
	except:
		print("Hydrogen bond between first and second hydration plot error")

os.system("mv *.txt *.jpg hydrogen-bond")

print("Hydrogen bond analysis done")
print("#------------------------------------------------------#")
print("\n")

# hydrogen bond lifetime in ions hydration shell

print("#------------------------------------------------------#")
print("Start hydrogen bond lifetime in ions hydration shell analysis")
os.system("mkdir hydrogenbond-lifetime")
os.system("cd hydrogenbond-lifetime")

for i in range(len(sel)):
	HBL_analysis=HBL(u,sel[i],sel[i],start,end,50)
	HBL_analysis.run()
	time=np.array([i for i in range(len(HBL_analysis.timeseries))]).\
	reshape(len(HBL_analysis.timeseries))
	HBLc=np.array([column[0] for column in HBL_analysis.timeseries]).\
	reshape(len(HBL_analysis.timeseries))
	HBLi=np.array([column[1] for column in HBL_analysis.timeseries]).\
	reshape(len(HBL_analysis.timeseries))
	data=np.zeros([len(HBL_analysis.timeseries),3])
	data[:,0]=time
	data[:,1]=HBLc
	data[:,2]=HBLi
	np.savetxt("HBL-%d-%d.txt"%(i,i),data,fmt="%.6f %.6f %.6f")
	print("%s HBL analysis done"%(sel[i]))
	
	try:
		plt.figure(1,figsize=(18,6))
		# HBL plot
		plt.subplot(121)
		plt.plot(time/10,HBLc,'g-',lw=2.0)
		plt.xlabel('Time (ns)')
		plt.ylabel("HBL Continuos")
		plt.plot("%s"%(sel[i]))
		plt.subplot(122)
		plt.plot(time/10,HBLi,'g-',lw=2.0)
		plt.xlabel('Time (ns)')
		plt.ylabel("HBL Intermitent")
		plt.plot("%s"%(sel[i]))
		plt.savefig("HBL-%d-%d.jpg"%(i,i),dpi=300)
		plt.clf()
	except:
		print("Hydrogen bond left time in first or second hydration shell plot error")

print("Hydrogen bond left time in first or second hydration shell analysis done")
print("#------------------------------------------------------#")
print("\n")
	
# hydrogen bond lifetime between ions first and second hydration shell

print("#------------------------------------------------------#")
print("Start hydrogen bond lifetime between ions first and second hydration shell analysis")

for i in range(0,len(sel),2):
	HBL_analysis=HBL(u,sel[i],sel[i+1],start,end,50)
	HBL_analysis.run()
	time=np.array([i for i in range(len(HBL_analysis.timeseries))]).\
	reshape(len(HBL_analysis.timeseries))
	HBLc=np.array([column[0] for column in HBL_analysis.timeseries]).\
	reshape(len(HBL_analysis.timeseries))
	HBLi=np.array([column[1] for column in HBL_analysis.timeseries]).\
	reshape(len(HBL_analysis.timeseries))
	data=np.zeros([len(HBL_analysis.timeseries),3])
	data[:,0]=time
	data[:,1]=HBLc
	data[:,2]=HBLi
	np.savetxt("HBL-%d-%d.txt"%(i,i+1),data,fmt="%.6f %.6f %.6f")
	print("%s and %s HBL analysis done"%(sel[i],sel[i+1]))
	
	# HBL plot
	try:
		plt.figure(1,figsize=(18,8))
		plt.subplot(121)
		plt.plot(time/10,HBLc,'g-',lw=2.0)
		plt.xlabel('Time (ns)')
		plt.ylabel("HBL Continuos")
		plt.title("%s\n %s\n"%(sel[i],sel[i+1]))
		plt.subplot(122)
		plt.plot(time/10,HBLi,'g-',lw=2.0)
		plt.xlabel('Time (ns)')
		plt.ylabel("HBL Intermitent")
		plt.title("%s\n %s\n"%(sel[i],sel[i+1]))
		plt.savefig("HBL-%d-%d.jpg"%(i,i+1),dpi=300)
		plt.clf()
	except:
		print("Hydrogen bond left time between first and second hydration shell plot error")
		
os.system("mv *.txt *.jpg hydrogenbond-lifetime")

print("Hydrogen bond lifetime in ions hydration shell analysis done")
print("#------------------------------------------------------#")
print("\n")
	
# Water Orientation Relaxation analysis

print("#------------------------------------------------------#")
print("Start Water Orientation Relaxation analysis")
os.system("mkdir Water-Orientation")
os.system("cd Water-Orientation")

for i in range(len(sel)):
	WOR_analysis=WOR(u,sel[i],start,end,50)
	WOR_analysis.run()
	data=np.zeros([len(WOR_analysis.timeseries),4])
	time=np.array([i for i in range(len(WOR_analysis.timeseries))]).\
	reshape(len(WOR_analysis.timeseries))
	data[:,0]=time
	data[:,1]=np.array([column[0] for column in WOR_analysis.timeseries]).\
	reshape(len(WOR_analysis.timeseries))
	data[:,2]=np.array([column[1] for column in WOR_analysis.timeseries]).\
	reshape(len(WOR_analysis.timeseries))
	data[:,3]=np.array([column[2] for column in WOR_analysis.timeseries]).\
	reshape(len(WOR_analysis.timeseries))
	
	np.savetxt("water-Orientation-%d-%d.txt"%(i,i),data,fmt="%.6f %.6f %.6f %.6f")
	print("%s water Orientation Relaxation analysis done"%(sel[i]))
	
	try:
		plt.figure(1,figsize=(18, 6))
		plt.subplot(131)
		plt.plot(data[:,0]/10,data[:,1],'g-',lw=2.0)
		plt.xlabel("Time")
		plt.ylabel("WOR")
		plt.title("WOR OH")
		
		plt.subplot(132)
		plt.plot(data[:,0]/10,data[:,2],'g-',lw=2.0)
		plt.xlabel("Time")
		plt.ylabel("WOR")
		plt.title("WOR HH")
		
		plt.subplot(133)
		plt.plot(data[:,0]/10,data[:,3],'g-',lw=2.0)
		plt.xlabel("Time")
		plt.ylabel("WOR")
		plt.title("WOR dip")
		
		plt.savefig("water-Orientation-%d-%d.jpg"%(i,i),dpi=300)
		plt.clf()
	except:
		print("Water Orientation Relaxation plot error")
		
os.system("mv *.txt *.jpg Water-Orientation")

print("Water Orientation Relaxation analysis done")
print("#------------------------------------------------------#")
print("\n")
	
# Water Angular Distribution analysis
# OH vector, HH vector,dipole vector
# P(cos(theta)) vs cos(theta)

print("#------------------------------------------------------#")
print("Start Water Angular Distribution analysis")
os.system("mkdir Water-Angular")
os.system("cd Water-Angular")

bins=30

for i in range(len(sel)):
	try:
		AD_analysis=AD(u,sel[i],bins,axis='z')
		AD_analysis.run()
		data=np.zeros([len(AD_analysis.graph[0][0]),4])
		data[:,0]=np.array([float(column.split()[0]) for column in AD_analysis.graph[0][:-1]]).\
		reshape(len(AD_analysis.graph[0][0]))
		data[:,1]=np.array([float(column.split()[1]) for column in AD_analysis.graph[0][:-1]]).\
		reshape(len(AD_analysis.graph[0][0]))
		data[:,2]=np.array([float(column.split()[1]) for column in AD_analysis.graph[1][:-1]]).\
		reshape(len(AD_analysis.graph[0][0]))
		data[:,3]=np.array([float(column.split()[1]) for column in AD_analysis.graph[2][:-1]]).\
		reshape(len(AD_analysis.graph[0][0]))
		np.savetxt("Water-Angular-Distribution-%d-%d.txt"%(i,i),data,fmt="%.6f %.6f %.6f %.6f")
		print("%s Water Angular Distribution analysis done"%(sel[i]))
	except:
		print("%s analysis error"%(sel[i]))	
	# AD OH
	try:
		plt.figure(1,figsize=(18, 6))
		plt.subplot(131)
		plt.plot(data[:,0],data[:,1],"g-",lw=2.0)
		plt.xlabel('cos theta')
		plt.ylabel('P(cos theta)')
		plt.title('PDF cos theta for OH')
		
		plt.subplot(132)
		plt.plot(data[:,0],data[:,2],"g-",lw=2.0)
		plt.xlabel('cos theta')
		plt.ylabel('P(cos theta)')
		plt.title('PDF cos theta for HH')
		
		plt.subplot(133)
		plt.plot(data[:,0],data[:,3],"g-",lw=2.0)
		plt.xlabel('cos theta')
		plt.ylabel('P(cos theta)')
		plt.title('PDF cos theta for dipole')
		plt.savefig("Water-Angular-Distribution-%d-%d.jpg"%(i,i),dpi=300)
		plt.clf()
	except:
		print("Water Angular Distribution plot error")
		
os.system("mv *.txt *.jpg Water-Angular")

print("Water Angular Distribution analysis done")
print("#------------------------------------------------------#")
print("\n")
	
# MSD analysis

print("#------------------------------------------------------#")
print("Start MSD analysis")
os.system("mkdir msd")
os.system("cd msd")

msdsel=[]
msdsel.append("name OW")
msdsel.append("resname %s"%(cation[0]))
msdsel.append("resname %s"%(anion[0]))

owcoord=np.zeros([u.trajectory.n_frames,u.select_atoms("%s"%(msdsel[0])).n_atoms,3])
cationcoord=np.zeros([u.trajectory.n_frames,u.select_atoms("%s"%(msdsel[1])).n_atoms,3])
anioncoord=np.zeros([u.trajectory.n_frames,u.select_atoms("%s"%(msdsel[2])).n_atoms,3])

for ts in u.trajectory:
	owcoord[ts.frame,:,:]=u.select_atoms("%s"%(msdsel[0])).positions
	cationcoord[ts.frame,:,:]=u.select_atoms("%s"%(msdsel[1])).positions
	anioncoord[ts.frame,:,:]=u.select_atoms("%s"%(msdsel[2])).positions
	
owmsd=np.zeros([u.trajectory.n_frames,4])
cationmsd=np.zeros([u.trajectory.n_frames,4])
anionmsd=np.zeros([u.trajectory.n_frames,4])
	
for i in range(u.select_atoms("%s"%(msdsel[0])).n_atoms):
	for frame in range(u.trajectory.n_frames):
		nt=u.trajectory.n_frames-frame
		for j in range(3):
			owmsd[frame,j] += np.sum((owcoord[1:nt,i,:]-owcoord[(1+frame):(nt+frame),i,:])**2)/float(nt)
		owmsd[frame,3] += sum(owmsd[frame,:3])
		
print("%s msd analysis done"%(msdsel[0]))

for i in range(u.select_atoms("%s"%(msdsel[1])).n_atoms):
	for frame in range(u.trajectory.n_frames):
		nt=u.trajectory.n_frames-frame
		for j in range(3):
			cationmsd[frame,j] += np.sum((cationcoord[1:nt,i,:]-cationcoord[(1+frame):(nt+frame),i,:])**2)/float(nt)
		cationmsd[frame,3] += sum(cationmsd[frame,:3])
		
print("%s msd analysis done"%(msdsel[1]))

for i in range(u.select_atoms("%s"%(msdsel[2])).n_atoms):
	for frame in range(u.trajectory.n_frames):
		nt=u.trajectory.n_frames-frame
		for j in range(3):
			anionmsd[frame,j] += np.sum((anioncoord[1:nt,i,:]-anioncoord[(1+frame):(nt+frame),i,:])**2)/float(nt)
		anionmsd[frame,3] += sum(anionmsd[frame,:3])
		
print("%s msd analysis done"%(msdsel[2]))
		
owmsd=owmsd/u.select_atoms("%s"%(msdsel[0])).n_atoms/100          # nm2
cationmsd=cationmsd/u.select_atoms("%s"%(msdsel[1])).n_atoms/100
anionmsd=anionmsd/u.select_atoms("%s"%(msdsel[2])).n_atoms/100

np.savetxt("msd-ow.txt",owmsd,fmt="%.6f %.6f %.6f %.6f")
np.savetxt("msd-cation.txt",cationmsd,fmt="%.6f %.6f %.6f %.6f")
np.savetxt("msd-anion.txt",anionmsd,fmt="%.6f %.6f %.6f %.6f")

try:
	plt.figure(1,figsize=(18,8))
	plt.subplot(131)
	labels=["x","y","z","all"]
	colors=["chartreuse","deepskyblue","b","r"]
	for i in range(4):
		plt.plot(owmsd[:,0],owmsd[:,i],color=colors[i],lw=2.0,label=labels[i])
	plt.xlabel('time')
	plt.ylabel('MSD')
	plt.legend(loc="best")
	plt.title('MSD for %s'%(msdsel[0]))
	
	plt.subplot(132)
	for i in range(4):
		plt.plot(cationmsd[:,0],cationmsd[:,i],color=colors[i],lw=2.0,label=labels[i])
	plt.xlabel('time')
	plt.ylabel('MSD')
	plt.legend(loc="best")
	plt.title('MSD for %s'%(msdsel[1]))
	
	plt.subplot(133)
	for i in range(4):
		plt.plot(anionmsd[:,0],anionmsd[:,i],color=colors[i],lw=2.0,label=labels[i])
	plt.xlabel('time')
	plt.ylabel('MSD')
	plt.legend(loc="best")
	plt.title('MSD for %s'%(msdsel[2]))
	
	plt.savefig("msd.jpg",dpi=300)
	plt.clf()

except:
	print("msd plot error")
		
os.system("mv *.txt *.jpg msd")

print("MSD analysis done")
print("#------------------------------------------------------#")
print("\n")
	
# Survival Probability == residence time

print("#------------------------------------------------------#")
print("Start Survival Probability analysis")
os.system("mkdir Survival-Probability")

for i in range(len(sel)):
	SP_analysis = SP(u,sel[i],start,end,20)
	SP_analysis.run()
	data=np.zeros([len(SP_analysis.timeseries),2])
	data[:,0]=np.array([i for i in range(len(SP_analysis.timeseries))]).\
	reshape(len(SP_analysis.timeseries))
	data[:,1]=np.array(SP_analysis.timeseries).reshape(len(SP_analysis.timeseries))
	np.savetxt("sp-%d.txt"%(i),data,fmt="%.6f %.6f")
	print("%s Survival Probability analysis done"%(sel[i]))
	
	try:
		plt.figure(1,figsize=(18,8))
		plt.plot(data[:,0],data[:,1],'g-',lw=2.0)
		plt.ylabel('SP')
		plt.title('Survival Probability of %s'%(sel[i]))
		plt.savefig("sp-%d.jpg"%(i),dpi=300)
		plt.clf()
	except:
		print("Survival Probability plot error")
		
os.system("mv *.txt *.jpg Survival-Probability")

print("Survival Probability analysis done")
print("#------------------------------------------------------#")
print("\n")
