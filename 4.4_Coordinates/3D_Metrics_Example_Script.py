#import what we need
import os; import chimera; from chimera import runCommand as rc; from chimera import replyobj; from chimera import Point; from chimera import selection; import Midas; from chimera import numpyArrayFromAtoms; import StructMeasure; from StructMeasure.Axes import axisManager; import csv; import numpy as np

# change to folder with data files
os.chdir("C:/Users/folder_where_the_PDBS_are")

#initialize table
table=[['fn','z_distance','z_short','Tilt','Translation_angle','Translation_distance','Rotation_angle','angleto30s']]
with open('chimera_attributes.csv', 'a') as csvfile:
		writer = csv.writer(csvfile)
		[writer.writerow(r) for r in table]
		
	
# gather the names of .pdb files in the folder
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

# loop through the files, opening, processing, and closing each in turn
for fn in file_names:
	replyobj.status("Processing " + fn) # show what file we're working on
	rc("open " + fn)
	rc("delete markers")
	#make receptor plane #select the top of each transmembrane# measure the center of the pocket (#1)#make a plane that mimics the mebrane (p1)
	rc("select #0:32,87,100,166,188,260,277"); rc("measure center sel mark true radius .3"); rc("define plane sel"); atoms=selection.currentAtoms(); numpyCorrds=numpyArrayFromAtoms(atoms); plane=StructMeasure.plane(numpyCorrds)
	#
	#select chemokine #measure axis of mass for the chemokine (a1)#measure center of mass for the chemokine (#2)#centriod for chemokine (c1)
	rc("select :.b"); rc("define axis color red  radius 2 sel"); rc("measure center sel mark true radius .3"); rc("define centroid mass true color blue radius 10 sel")
	#
	#
	z_distance=Midas.distance('#1:1' '#2:1')#distance center of mass of chemokine to center of pocket at membrane
	z_short=Midas.distance('#1',['c1','p1'])#distance from center of mass of chemokine to membrane 
	Tilt=Midas.angle('#1',['a1','p1'])#angle of tilt from membrane to chemokine axis
	#
	## select C in Nterm of receptor as 0 for polar coordinates (#1:2)#convert to atoms#convert to Point#project Point onto plane
	rc("select #0:20@CA"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom)
	#
	#delete point point so all future is the same
	rc("select #1"); rc("delete sel")
	#
	#(1:2)select the top of each transmembrane# measure the center of the pocket (#1)
	rc("select #0:32,87,100,166,188,260,277"); rc("measure center sel mark true radius .3"); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")
	#
	# select center of mass of chemokine (#1:3)#convert to atoms#convert to Point#project Point onto plane#project Point onto plane
	rc("select #2"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")
	#
	# select F at base of ecl2 (#1:4)#convert to atoms
	rc("select #0:166@CA"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")
	#
	#select CYS at end of ecl1 (#1:5)
	rc("select #0:101@CA"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")
	#
	#
	rc("select #1:1,4"); rc("define axis color red  radius 1 sel");
	Translation_angle=Midas.angle('#1:3 #1:1 #1:4')
	rc("define axis p1")
	#
	##1:6 = receptor_90_guideline
	axes = axisManager.axes; axes.sort(lambda a2, a3: cmp(a2.number, a3.number)); ax3=axes[2]; cx, cy, cz = ax3.center; dx, dy, dz = ax3.direction; vec3=[dx,dy,dz]; d0x, d0y, d0z = axes[1].direction; vec2=[d0x,d0y,d0z]; ax4=np.cross(vec2, vec3)
	if ax4[0] < 0: ax4 = ax4*-1
	rc("cofr " + str((ax4[0]*10)+cx) + "," + str((ax4[1]*10)+cy) + "," + str((ax4[2]*10)+cz) + " coord #0; ac mc")##1:8
	Translation_angle_90=Midas.angle('#1:6 #1:1 #1:3')
	if Translation_angle_90 < 90: Translation_angle = 360-Translation_angle
	Translation_distance=Midas.distance('#1:1' '#1:3')
	#
	#p2
	rc("select :.b"); rc("define plane sel"); rc("define axis color red  radius 2 p2")#a2
	axes = axisManager.axes; axes.sort(lambda a1, a4: cmp(a1.number, a4.number)); ax4=axes[3]
	#end1= [a1.direction * ext + a1.center for ext in a1.extents]
	#cx, cy, cz = Point(end1)
	cx, cy, cz = ax4.center; dx, dy, dz = ax4.direction; vec1=[dx,dy,dz]
	#
	d0x, d0y, d0z = axes[0].direction; vec0=[d0x,d0y,d0z]; ax2=np.cross(vec1, vec0)
	#
	rc("cofr " + str(dx+cx) + "," + str(dy+cy) + "," + str(dz+cz) + " coord #0; ac mc") ##1:7
	rc("cofr " + str((dx*-1)+cx) + "," + str((dy*-1)+cy) + "," + str((dz*-1)+cz) + " coord #0; ac mc")##1:8
	rc("cofr " + str((ax2[0])+cx) + "," + str((ax2[1])+cy) + "," + str((ax2[2])+cz) + " coord #0; ac mc")##1:9
	rc("cofr " + str((ax2[0]*-1)+cx) + "," + str((ax2[1]*-1)+cy) + "," + str((ax2[2]*-1)+cz) + " coord #0; ac mc")##1:10
	rc("select #1:7-10"); rc("define plane radius 20 sel")#p3
	#
	atoms=selection.currentAtoms(); numpyCorrds=numpyArrayFromAtoms(atoms); plane=StructMeasure.plane(numpyCorrds)
	#
	#F before ECL2 as 0point reference #1:11
	rc("select #0:166@CA"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")#project Point onto plane
	#
	#select 30s loop ##1:12
	rc("select #0:350@CA"); x=selection.currentAtoms();	myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")#project Point onto plane
	#
	#select cys in chemokien B3 ##1:13
	rc("select #0:367@CA"); x=selection.currentAtoms(); myatom=x[0].coord(); nearest=plane.nearest(myatom); rc("cofr " + str(nearest[0]) + "," + str(nearest[1]) + "," + str(nearest[2]) + " coord #0; ac mc")#project Point onto plane
	#
	Rotation_angle=Midas.angle('#1:11 #2 #1:13') #angle of chemokoine roatation
	angleto30s=Midas.angle('#1:11 #2 #1:12') #angle of chemokoine roatation to 30s loop
	rc("select #1:11 #2"); rc("define axis color red  radius 1 sel")
	rc("define axis p3")
	axes = axisManager.axes; axes.sort(lambda a5, a6: cmp(a5.number, a6.number)); ax6=axes[5]; cx, cy, cz = ax6.center; dx, dy, dz = ax6.direction; vec3=[dx,dy,dz]; d0x, d0y, d0z = axes[4].direction; vec2=[d0x,d0y,d0z]; ax4=np.cross(vec2, vec3)
	if ax4[0] < 0: ax4 = ax4*-1
	rc("cofr " + str((ax4[0]*10)+cx) + "," + str((ax4[1]*10)+cy) + "," + str((ax4[2]*10)+cz) + " coord #0; ac mc")##1:8
	Rotation_angle_90=Midas.angle('#1:14 #2 #1:13')
	if Rotation_angle_90 < 90: Rotation_angle = 360-Rotation_angle
	angleto30s_90=Midas.angle('#1:14 #2 #1:12')
	if angleto30s_90 < 90: angleto30s = 360-angleto30s
	#
	table=[[fn,z_distance,z_short,Tilt,Translation_angle,Translation_distance,Rotation_angle,angleto30s]]
	with open('chimera_attributes.csv', 'a') as csvfile:
		writer = csv.writer(csvfile)
		[writer.writerow(r) for r in table]
	rc("close session")
rc("stop now")

