#import os
#import sys
import decimal
import numpy as np
import math as mt
import pandas as pd



print("|-----------------------------------------------------| \r\n|----------------------B.E.L.L.O----------------------| \r\n|---------Bond Element Lattice Locality Order---------|\r\n|-----------------------------------------------------|")

f = pd.read_fwf('traj.xyz',header=None)
f= f.fillna("x") #fill the empty spaces with X so that we can clean/organize them before analysis
file = np.array(f)

Natom=int(f[0][0])

trhld=3.0#float(input('Input atomic distance threshold (Ang): '))

tlrnc=0.10#float(input('Input atomic distance tolerance (Ang): '))

#celldm=celldmx=celldmy=celldmz=20.644848284633
#celldm=celldmx=celldmy=celldmz=20.644848284633#float(input('If your unitcell is cubic, then enter the lattice constant in Ang: '))
celldmx=20.644848284633	
celldmy=20.644848284633
celldmz=20.644848284633

lfole=len(file)
print("File lenght is: ",lfole)

#------------------------------------
#-----variables and empy lists-------
#------------------------------------
fa= np.empty((0,4), int)
zerocenter=[0,0,0]
finalq=[]
atomnumbers=[] #used for finding atom numbers for pdos summation
coordinates=[] #human readable coordinates
xyzcoordinates= [] # final xyz coordinates
pdbcoords=[]# Pdb format output
angle=[] # angle distribution
localstats=[] # output for statistics on local orders
localstats.append(["3-FOLD","4-FOLD","TETRAHEDRAL","5-FOLD","OCTAHEDRAL","Total"])
t1PDB= "PDB coordinates"
t2PDB="-------BELLO-------"
t3PDB="CRYST1   %f   %f   %f  90.00  90.00  90.00 P1          1" % (celldmx, celldmy, celldmz)
pdbcoords.append(t1PDB)
pdbcoords.append(t2PDB)
pdbcoords.append(t3PDB)
qframe=[]
q3flframe=[]
q4flframe=[]
qtetframe=[]
q5flframe=[]
q6flframe=[]
tetatot=[]

#------------------------------------
#------boundry condition #1----------
#------------------------------------

def distance(a, b):
	dx = abs(a[0] - b[0])
	x = min(dx, abs(A - dx))

	dy = abs(a[1] - b[1])
	y = min(dy, abs(B - dy))

	dz = abs(a[2] - b[2])
	z = min(dz, abs(C - dz))

	return mt.sqrt(x**2 + y**2 + z**2)

A=celldmx
B=celldmy
C=celldmz

#------------------------------------
#-----------Dot function-------------
#------------------------------------

def dott(a,b,c):
#a and b are the bonds
#c is the center atom
	celld=[celldmx,celldmy,celldmz]
	tlist=[a,b]
	#np.linalg.norm(j-c))
	for j in tlist:
		d=(mt.sqrt((j[0]-c[0])**2 + (j[1]-c[1])**2 + (j[2]-c[2])**2))
		if ( d != distance(j,c)):
			for i in range(0,3):
				if abs(j[i]-c[i])>(celld[i]/2):
					if c[i] < j[i]:
						j[i]=(j[i] - celld[i])
					elif c[i] > j[i]:
						j[i]=(j[i] + celld[i])

	return sum((tlist[0]-c)*(tlist[1]-c))

#------------------------------------
#----------Duplicate handler---------
#------------------------------------

def rep(a):
	binn=[]
	binn.append(a[0])
	for i in range(1,len(a)):
		if np.any([a[i][1:4] != x[1:4] for x in binn]):
			binn.append(a[i])

	return binn

#------------------------------------
#------boundry condition #2----------
#------------------------------------
for i in range(lfole):
	if (file[i,1] != 'x'):
		if file[i,1]>celldmx or file[i,1]<0:
			file[i,1]= file[i,1]%celldmx
		if file[i,2]>celldmy or file[i,2]<0:
			file[i,2]= file[i,2]%celldmy
		if file[i,3]>celldmz or file[i,3]<0:
			file[i,3]= file[i,3]%celldmz
		fa= np.vstack((fa, file[i]))

print("Boundry Condition is done! \r\nCalculating local orders:")
#----------------------------------
for N in range(0,len(fa),Natom):
#-variables for each frame-----
	fo=[]
	fi=[]
	fo=np.copy(fa[N:N+Natom])
	lfo=len(fo)
	q=[]
	q3fl=[]
	q4fl=[]
	qtet=[]
	q5fl=[]
	q6fl=[]
	aq=[]
	mcoordinates= [] # temporary xyz coordinates
	nlocal=0 # number of local orders
	n3fl=0   # number of 3-FOLDs
	n4fl=0   # number of 4-FOLDs
	ntet=0   # number of TETRAHEDRALs
	n5fl=0   # number of 5-FOLDs
	noct=0   # number of OCTAHEDRALs
	condition='false' #this is an on-off switch for when we calculate the q, to tell the code to jump back up to selecting b1
	for c in range (0,lfo):  #in here we select the point which will be the center of tetrahedral
		templist=[]
		condition='false'
		Gecounter=0 # this one counts the number of Ge atoms so in the end we can check to see whether they are in the center or not
		nbond=0 # number of bonds
		if condition=='false':
			for b1 in range(0,lfo): #selects fisrt bond (b1)
				if c!=b1 and condition=='false':   #checks not to select the same atom as center atom
					a=np.copy(fo[c,1:4])
					b=np.copy(fo[b1,1:4])
					ab1=np.copy(fo[b1,1:4])
					db1=distance(a,b)
					if (trhld-tlrnc)<db1<(trhld+tlrnc):
						nbond += 1 # now we have 1 bond
						if condition=='false':
							for b2 in range(0,lfo):   #selects second bond and checks to see if it is different from center and b1
								if c!=b2 and b1!=b2 and condition=='false':
									b=np.copy(fo[b2,1:4])
									ab2=np.copy(fo[b2,1:4])
									db2=distance(a,b)
									if abs(db1-db2) < tlrnc:  #if the bond lenght of both b1 and b2 are the same, then calculate the COS
										nbond += 1 # now we have 2 bonds

										dot = dott(ab1,ab2,a)
										cosb12 = dot / (db1 * db2)
										teta = mt.degrees(mt.acos(cosb12))

										if condition=='false':
											for b3 in range(0,lfo):
												if b3!=c and b3!=b2 and b3!=b1 and condition=='false':
													b=np.copy(fo[b3,1:4])
													ab3=np.copy(fo[b3,1:4])
													db3=distance(a,b)
													if abs(db3-db1)<=tlrnc:
														nbond += 1
														
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b2,0])) #teta from previous bond but should be appended here, otherwise we will have teta for 2folds
														dot = dott(ab3,ab2,a)
														cosb23 = dot / (db3 * db2)
														teta = mt.degrees(mt.acos(cosb23))
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b3,0]))

														dot = dott(ab1,ab3,a)
														cosb13 = dot / (db3 * db1)
														teta = mt.degrees(mt.acos(cosb13))
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b3,0]))

														if condition=='false':
															for b4 in range(0,lfo):
																if c!=b4 and b1!=b4 and b2!=b4 and b3!=b4 and condition=='false':
																	b=np.copy(fo[b4,1:4])
																	ab4=np.copy(fo[b4,1:4])
																	db4=distance(a,b)
																	if abs(db4-db1)<=tlrnc:
																		nbond += 1
																		dot = dott(ab1,ab4,a)
																		cosb14 = dot / (db1 * db4)
																		teta = mt.degrees(mt.acos(cosb14))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b4,0]))

																		dot = dott(ab2,ab4,a)
																		cosb24 = dot / (db2 * db4)
																		teta = mt.degrees(mt.acos(cosb24))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b4,0]))

																		dot = dott(ab4,ab3,a)
																		cosb34 = dot / (db3 * db4)
																		teta = mt.degrees(mt.acos(cosb34))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b4,0]))

																		if condition=='false':
																			for b5 in range(0,lfo):
																				if c!=b5 and b1!=b5 and b2!=b5 and b3!=b5 and b4!=b5 and condition=='false':
																					b=np.copy(fo[b5,1:4])
																					ab5=np.copy(fo[b5,1:4])
																					db5=distance(a,b)
																					if abs(db5-db1)<=tlrnc:
																						nbond += 1

																						dot = dott(ab1,ab5,a)
																						cosb15 = dot / (db1 * db5)
																						teta = mt.degrees(mt.acos(cosb15))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b5,0]))

																						dot = dott(ab2,ab5,a)
																						cosb25 = dot / (db2 * db5)
																						teta = mt.degrees(mt.acos(cosb25))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b5,0]))

																						dot = dott(ab3,ab5,a)
																						cosb35 = dot / (db3 * db5)
																						teta = mt.degrees(mt.acos(cosb35))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b5,0]))

																						dot = dott(ab4,ab5,a)
																						cosb45 = dot / (db4 * db5)
																						teta = mt.degrees(mt.acos(cosb45))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b4,0],fo[c,0],fo[b5,0]))

																						if condition=='false':
																							for b6 in range(0,lfo):
																								if c!=b6 and b1!=b6 and b2!=b6 and b3!=b6 and b4!=b6 and b5!=b6 and condition=='false':
																									b=np.copy(fo[b6,1:4])
																									ab6=np.copy(fo[b6,1:4])
																									db6=distance(a,b)
																									if abs(db6-db1)<=tlrnc:
																										nbond += 1

																										dot = dott(ab1,ab6,a)
																										cosb16 = dot / (db1 * db6)
																										teta = mt.degrees(mt.acos(cosb16))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b6,0]))

																										dot = dott(ab2,ab6,a)
																										cosb26 = dot / (db2 * db6)
																										teta = mt.degrees(mt.acos(cosb26))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b6,0]))

																										dot = dott(ab3,ab6,a)
																										cosb36 = dot / (db3 * db6)
																										teta = mt.degrees(mt.acos(cosb36))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b6,0]))

																										dot = dott(ab4,ab6,a)
																										cosb46 = dot / (db4 * db6)
																										teta = mt.degrees(mt.acos(cosb46))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b4,0],fo[c,0],fo[b6,0]))

																										dot = dott(ab5,ab6,a)
																										cosb56 = dot / (db5 * db6)
																										teta = mt.degrees(mt.acos(cosb56))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b5,0],fo[c,0],fo[b6,0]))
																										condition='true'
																									else:
																										continue
																								if condition=='true':
																									break
																							condition=='true'
																							break
																						else:
																							continue
																					else:
																						continue
																				if condition=='true':
																					break
																			condition='true'
																			break
																		else:
																			continue
																	else:
																		continue
																if condition=='true':
																	break
															condition='true'
															break
													else:
														continue
												if condition=='true':
													break
											condition='true'
											break
										else:
											continue
									else:
										continue
								if condition=='true':
									break
							condition='true'
							break
						else:
							continue
					else:
						continue
				if condition=='true':
					break
			else:
				condition='true'

			if nbond==3:
				nlocal += 1
				n3fl += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb23 + 1/3)**2))
				q3fl.append(q)
				aq.append(q) # aq is the list of q values for each frame
				templist= [c,b1,b2,b3]
				coordinates.append(["3-FOLD","-","q is:",q])
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('3-FOLD',c+1,b1+1,b2+1,b3+1,0,0,0)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4]) 
					coordinates.append(fo[i,0:4])
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'3FL','A',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)
				coordinates.append([" "," "," "," "])
				condition='true' # if the code reaches here, we want it to jump back to selecting b1
			if nbond==4:
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2))
				aq.append(q)
				templist= [c,b1,b2,b3,b4]
				if 1>= q >=0.85:
					qtet.append(q)
					nlocal += 1
					ntet += 1
					templist= [c,b1,b2,b3,b4]
					coordinates.append(["TETRAHEDRAL","-","q is:",q])
					atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('TETHDL',c+1,b1+1,b2+1,b3+1,b4+1,0,0)
					atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
					for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
						mcoordinates.append(fo[i,0:4])
						coordinates.append(fo[i,0:4])
						pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'TET','B',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" "," "," "," "])
				elif   q <=0.85:
					q4fl.append(q)
					nlocal += 1
					n4fl += 1
					coordinates.append(["4-FOLD","-","q is:",q])
					atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('4-FOLD',c+1,b1+1,b2+1,b3+1,b4+1,0,0)
					atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
					for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
						mcoordinates.append(fo[i,0:4])
						coordinates.append(fo[i,0:4])
						pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'4FL','C',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" "," "," "," "])
				condition='true'
			if nbond==5:
				nlocal += 1
				n5fl += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2 + (cosb15 + 1/3)**2 + (cosb25 + 1/3)**2 + (cosb35 + 1/3)**2 + (cosb45 + 1/3)**2))
				q5fl.append(q)
				aq.append(q)
				templist= [c,b1,b2,b3,b4,b5]
				coordinates.append(["5-FOLD","-","q is:",q])
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('5-FOLD',c+1,b1+1,b2+1,b3+1,b4+1,b5+1,0)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4])
					coordinates.append(fo[i,0:4])
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'5FL','D',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)
				coordinates.append([" "," "," "," "])
				condition='true'
			if nbond==6:
				nlocal += 1
				noct += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2 + (cosb15 + 1/3)**2 + (cosb25 + 1/3)**2 + (cosb35 + 1/3)**2 + (cosb45 + 1/3)**2 + (cosb16 + 1/3)**2 + (cosb26 + 1/3)**2 + (cosb36 + 1/3)**2 + (cosb46 + 1/3)**2 + (cosb56 + 1/3)**2))
				q6fl.append(q)
				aq.append(q)
				templist= [c,b1,b2,b3,b4,b5,b6]
				coordinates.append(["OCTAHEDRAL","-","q is:",q])
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('OCTHDL',c+1,b1+1,b2+1,b3+1,b4+1,b5+1,b6+1)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4])
					coordinates.append(fo[i,0:4])
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'OCT','E',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)
				coordinates.append([" "," "," "," "])
				condition='true'
		else:
			continue
	Endf= "End of frame %d " % (N/Natom)
	for x in rep(mcoordinates):
		xyzcoordinates.append(x) # xyz coords of local orders
	xyzcoordinates.append(["\r\n",Endf,"\r\n"," "])
	atomnumbers.append(Endf)
	stat= " \r\n Number of 3-FOLD, 4-FOLD, TETRAHEDRAL, 5-FOLD, OCTAHEDRAL: %d, %d, %d, %d, %d" % (n3fl, n4fl, ntet, n5fl, noct)
	locq= "\r\n Number of local orders are: %d " % (len(aq))
	coordinates.append([Endf,stat,locq,"\r\n"])
	localstats.append([n3fl,n4fl,ntet,n5fl,noct,nlocal])
	pdbcoords.append(Endf)
	aq=[round(num,2) for num in aq]
	qframe.append(aq)
	q3fl=[round(num,2) for num in q3fl]
	q3flframe.append(q3fl)
	q4fl=[round(num,2) for num in q4fl]
	q4flframe.append(q4fl)
	qtet=[round(num,2) for num in qtet]
	qtetframe.append(qtet)
	q5fl=[round(num,2) for num in q5fl]
	q5flframe.append(q5fl)
	q6fl=[round(num,2) for num in q6fl]
	q6flframe.append(q6fl)
	print(Endf)


np.savetxt('output-q3fl.txt', q3flframe, delimiter=',', fmt="%s")
np.savetxt('output-q4fl.txt', q4flframe, delimiter=',', fmt="%s")
np.savetxt('output-qtet.txt', qtetframe, delimiter=',', fmt="%s")
np.savetxt('output-q5fl.txt', q5flframe, delimiter=',', fmt="%s")
np.savetxt('output-q6fl.txt', q6flframe, delimiter=',', fmt="%s")
np.savetxt('output-q-total.txt', qframe, delimiter=',', fmt="%s")
np.savetxt('output-human-readable-coords.txt', coordinates, delimiter=' ', fmt="%s")
np.savetxt('output-atom-number.txt',atomnumbers,delimiter=',',fmt="%s")
np.savetxt('output-xyz-coords.txt', xyzcoordinates, delimiter=' ', fmt="%s")
np.savetxt('output-pdb-coords-m.txt', pdbcoords, fmt="%s")
np.savetxt('output-local-statistics.txt', localstats, fmt="%s")
np.savetxt('output-angle-distribution.txt', tetatot, fmt="%s")
print("Done!")

