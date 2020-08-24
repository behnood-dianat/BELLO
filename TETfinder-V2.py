#import os
#import sys
import decimal
import numpy as np
import math as mt
import pandas as pd



print("|-----------------------------------------------------| \r\n|----------------------B.E.L.L.O----------------------| \r\n|---------Bond Element Lattice Locality Order---------|\r\n|-----------------------------------------------------|")
#Natom= int(input('Enter the number of atoms: '))
Natom=336
f = pd.read_fwf('traj.xyz',header=None)
f= f.fillna("x") #fill the empty spaces with X so that we can clean/organize them before analysis
file = np.array(f)

#emptylist=[" "," "," "," "]
#celldm=20.644848284633
#celldm=float(input('If your unitcell is cubic, then enter the lattice constant in Ang: '))
celldmx=19.12116
celldmy=19.38795
celldmz=21.45236
lfile=len(file)
print("File lenght is: ",lfile)

#------------------------------------
#-----variables and empy lists-------
#------------------------------------
fa= np.empty((0,4), int)
zerocenter=[0,0,0]
finalq=[]
atomnumbers=[] #used for finding atom numbers for pdos summation
coordinates=[]#= np.empty((0,4), int) #human readable coordinates
mcoordinates= [] #np.empty([5000000,4], int) #standard xyz coordinates
pdbcoords=[]#= np.empty([], str) # Pdb format output
angle=[]#= np.empty([],str) # angle distribution
localstats=[] #np.empty([lfifile,6], int) # output for statistics on local orders
localstats.append(["3-FOLD","4-FOLD","TETRAHEDRAL","5-FOLD","OCTAHEDRAL","Total"])# np.append(localstats, [["3-FOLD","4-FOLD","TETRAHEDRAL","5-FOLD","OCTAHEDRAL","Total"]], axis=0)
tPDB= "PDB coordinates"
pdbcoords.append([tPDB])#= np.vstack((pdbcoords, tPDB))
qframe=[]
q3flframe=[]
q4flframe=[]
qtetframe=[]
q5flframe=[]
q6flframe=[]
tetatot=[]
#------------------------------------
#------boundry condition-------------
#------------------------------------
for i in range(lfile):
	if (file[i,1] != 'x'):
		if file[i,1]>celldmx or file[i,1]<0:
			file[i,1]= file[i,1]%celldmx
		if file[i,2]>celldmy or file[i,2]<0:
			file[i,2]= file[i,2]%celldmy
		if file[i,3]>celldmz or file[i,3]<0:
			file[i,3]= file[i,3]%celldmz
		fa= np.vstack((fa, file[i]))
#		print(file[i])
#		llf=len(fa)
#		print(fa[llf-1])
print("Boundry Condition is done! \r\nCalculating local orders:")
#----------------------------------
for N in range(0,len(fa),Natom):
#-variables for each frame-----
	fo=[]
	fi=[]
	fo=np.copy(fa[N:N+Natom])
	fi=np.copy(fa[N:N+Natom])
	lfo=len(fo)
	lfi=len(fi)
	q=[]
	q3fl=[]
	q4fl=[]
	qtet=[]
	q5fl=[]
	q6fl=[]
	aq=[]
	nlocal=0 # number of local orders
	n3fl=0   # number of 3-FOLDs
	n4fl=0   # number of 4-FOLDs
	ntet=0   # number of TETRAHEDRALs
	n5fl=0   # number of 5-FOLDs
	noct=0   # number of OCTAHEDRALs
	condition='false' #this is an on-off switch for when we calculate the q, to tell the code to jump back up to selecting b1
	for c in range (0,lfo):  #in here we select the point which will be the center of tetrahedral
		#print("center is", fo[c])
		#print("fo is:", fo)
		if (fo[c,1:4] != zerocenter).any(): #checks if the center is [0 0 0] or not, otherwise it will shift the whole coordinates so that selected atom is [0 0 0]
			for i in range(0,lfi):
				fi[i,1:4]=fo[i,1:4]-fo[c,1:4]
				#print("changed coordinate f%i " %(i), fi[i])
		else:
			pass

		center= fi[c,1:4]
		templist=[]
		#print("center is: ", c, templist)
		condition='false'
		Gecounter=0 # this one counts the number of Ge atoms so in the end we can check to see whether they are in the center or not
		nbond=0 # number of bonds
		if condition=='false':
			for b1 in range(0,lfi): #selects fisrt bond (b1)
				if c!=b1 and condition=='false':   #checks not to select the same atom as center atom
					#print("b1 is:", b1)
					ab1=fi[b1,1:4]          #put a line of array (xyz) into a list
					distance=(((center[0]-ab1[0])**2)+((center[1]-ab1[1])**2)+((center[2]-ab1[2])**2))**0.5
					db1=distance
					if db1<3.0: #if the bound lenght is bigger than 3.2, we dont want it
						nbond += 1 # now we have 1 bond
						if condition=='false':
							for b2 in range(0,lfi):   #selects second bond and checks to see if it is different from center and b1
								if c!=b2 and b1!=b2 and condition=='false':
									#print("b2 is:", b2)
									ab2=fi[b2,1:4]
									distance=(((center[0]-ab2[0])**2)+((center[1]-ab2[1])**2)+((center[2]-ab2[2])**2))**0.5
									db2=distance
									if abs(db1-db2)< 2e-1:  #if the bond lenght of both b1 and b2 are the same, then calculate the COS
										nbond += 1 # now we have 2 bonds
										#print("b2")

										dot = sum(ab1*ab2)
										#dot = np.dot(ab1, ab2)
										normb1 = db1
										normb2 = db2
										cosb12 = dot / (normb1 * normb2)
										teta = mt.degrees(mt.acos(cosb12))

										#print("cos 12 is:", cosb12) #cos12 means Cos angle-b1-b2
										if condition=='false':
											for b3 in range(0,lfi):
												if b3!=c and b3!=b2 and b3!=b1 and condition=='false':
													#print("b3 is:", b3)
													ab3=fi[b3,1:4]
													distance=(((center[0]-ab3[0])**2)+((center[1]-ab3[1])**2)+((center[2]-ab3[2])**2))**0.5
													db3=distance
													if abs(db3-db1)<2e-1:
														nbond += 1
														#print("b3")
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b2,0])) #teta from previous bond but should be appended here, otherwise we will have teta for 2folds
														dot = sum(ab3*ab2)
														#dot = np.dot(ab3, ab2)
														normb3 = db3
														normb2 = db2
														cosb23 = dot / (normb3 * normb2)
														teta = mt.degrees(mt.acos(cosb23))
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b3,0]))

														dot = sum(ab3*ab1)
														#dot = np.dot(ab3, ab1)
														normb3 = db3
														normb1 = db1
														cosb13 = dot / (normb3 * normb1)
														teta = mt.degrees(mt.acos(cosb13))
														tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b3,0]))

														#print("cos 23-13 is:", cosb23, cosb13)
														if condition=='false':
															for b4 in range(0,lfi):
																if c!=b4 and b1!=b4 and b2!=b4 and b3!=b4 and condition=='false':
																	#print("b4 is:", b4)
																	ab4=fi[b4,1:4]
																	distance=(((center[0]-ab4[0])**2)+((center[1]-ab4[1])**2)+((center[2]-ab4[2])**2))**0.5
																	db4=distance
																	if abs(db4-db1)<2e-1:
																		nbond += 1
																		#print("b4")
																		dot = sum(ab4*ab1)
																		#dot = np.dot(ab4, ab1)
																		normb4 = db4
																		normb1 = db1
																		cosb14 = dot / (normb1 * normb4)
																		teta = mt.degrees(mt.acos(cosb14))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b4,0]))

																		dot = sum(ab4*ab2)
																		#dot = np.dot(ab4, ab2)
																		normb4 = db4
																		normb2 = db2
																		cosb24 = dot / (normb2 * normb4)
																		teta = mt.degrees(mt.acos(cosb24))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b4,0]))

																		dot = sum(ab4*ab3)
																		#dot = np.dot(ab4, ab3)
																		normb3 = db3
																		normb4 = db4
																		cosb34 = dot / (normb3 * normb4)
																		teta = mt.degrees(mt.acos(cosb34))
																		tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b4,0]))

																		if condition=='false':
																			for b5 in range(0,lfi):
																				if c!=b5 and b1!=b5 and b2!=b5 and b3!=b5 and b4!=b5 and condition=='false':
																					#print("b5 is:", b5)
																					ab5=fi[b5,1:4]
																					distance=(((center[0]-ab5[0])**2)+((center[1]-ab5[1])**2)+((center[2]-ab5[2])**2))**0.5
																					db5=distance
																					if abs(db5-db1)<2e-1:
																						nbond += 1
																						#print("b5")

																						dot = sum(ab5*ab1)
																						normb5 = db5
																						normb1 = db1
																						cosb15 = dot / (normb1 * normb5)
																						teta = mt.degrees(mt.acos(cosb15))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b5,0]))

																						dot = sum(ab5*ab2)
																						normb5 = db5
																						normb2 = db2
																						cosb25 = dot / (normb2 * normb5)
																						teta = mt.degrees(mt.acos(cosb25))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b5,0]))

																						dot = sum(ab5*ab3)
																						normb5 = db5
																						normb3 = db3
																						cosb35 = dot / (normb3 * normb5)
																						teta = mt.degrees(mt.acos(cosb35))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b5,0]))

																						dot = sum(ab4*ab5)
																						normb5 = db5
																						normb4 = db4
																						cosb45 = dot / (normb4 * normb5)
																						teta = mt.degrees(mt.acos(cosb45))
																						tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b4,0],fo[c,0],fo[b5,0]))

																						if condition=='false':
																							for b6 in range(0,lfi):
																								#print("looping b6", b6)
																								if c!=b6 and b1!=b6 and b2!=b6 and b3!=b6 and b4!=b6 and b5!=b6 and condition=='false':
																									#print("b6 is:", b6)
																									ab6=fi[b6,1:4]
																									distance=(((center[0]-ab6[0])**2)+((center[1]-ab6[1])**2)+((center[2]-ab6[2])**2))**0.5
																									db6=distance
																									if abs(db6-db1)<2e-1:
																										nbond += 1
																										#print("b6")

																										dot = sum(ab6*ab1)
																										normb6 = db6
																										normb1 = db1
																										cosb16 = dot / (normb1 * normb6)
																										teta = mt.degrees(mt.acos(cosb16))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b1,0],fo[c,0],fo[b6,0]))

																										dot = sum(ab6*ab2)
																										normb6 = db6
																										normb2 = db2
																										cosb26 = dot / (normb2 * normb6)
																										teta = mt.degrees(mt.acos(cosb26))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b2,0],fo[c,0],fo[b6,0]))

																										dot = sum(ab6*ab3)
																										normb6 = db6
																										normb3 = db3
																										cosb36 = dot / (normb3 * normb6)
																										teta = mt.degrees(mt.acos(cosb36))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b3,0],fo[c,0],fo[b6,0]))

																										dot = sum(ab4*ab6)
																										normb6 = db6
																										normb4 = db4
																										cosb46 = dot / (normb4 * normb6)
																										teta = mt.degrees(mt.acos(cosb46))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b4,0],fo[c,0],fo[b6,0]))

																										dot = sum(ab5*ab6)
																										normb6 = db6
																										normb5 = db5
																										cosb56 = dot / (normb5 * normb6)
																										#print("dot db5 db6 cosb56",ab5,ab6,dot,db6,db5,cosb56)
																										#print(fi)
																										teta = mt.degrees(mt.acos(cosb56))
																										tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta,fo[b5,0],fo[c,0],fo[b6,0]))
																										#print(tetatot)
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
			#if nbond==3:
			#	templist= [c,b1,b2,b3]
			#	for i in templist:
			#		if fo[i,0]=='Ge':
			#			Gecounter += 1
			#if nbond==4:
			#	templist= [c,b1,b2,b3,b4]
			#	for i in templist:
			#		if fo[i,0]=='Ge':
			#			Gecounter += 1
			#if nbond==5:
			#	templist= [c,b1,b2,b3,b4,b5]
			#	for i in templist:
			#		if fo[i,0]=='Ge':
			#			Gecounter += 1
			#if nbond==6:
			#	templist= [c,b1,b2,b3,b4,b5,b6]
			#	for i in templist:
			#		if fo[i,0]=='Ge':
			#			Gecounter += 1
			#if (Gecounter >= 2 and fo[c,0]=='Ge') or (Gecounter <= 1):
			if nbond==3:
				nlocal += 1
				n3fl += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb23 + 1/3)**2))
				q3fl.append(q)
				aq.append(q) # aq is the array of q values for each frame (270atoms)
				templist= [c,b1,b2,b3]
				coordinates.append(["3-FOLD","-","q is:",q])#= np.append(coordinates, [["3-FOLD","-","-","-"]], axis=0)
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('3-FOLD',c+1,b1+1,b2+1,b3+1,0,0,0)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4]) #np.vstack((mcoordinates, fo[i,0:4]))
					coordinates.append(fo[i,0:4])#= np.vstack((coordinates, fo[i,0:4]))
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'3FL','A',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)#= np.vstack((pdbcoords, pdbtemp))
				coordinates.append([" "," "," "," "])#= np.append(coordinates, [[" "," "," "," "]], axis=0)
				condition='true' # if the code reaches here, we want it to jump back to selecting b1
			if nbond==4:
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2))
				#print(cosb12,cosb13,cosb14,cosb23,cosb24,cosb34)
				aq.append(q)
				templist= [c,b1,b2,b3,b4]
				#print("q n4:", q)
				if 1>= q >=0.85:
					qtet.append(q)
					nlocal += 1
					ntet += 1
					templist= [c,b1,b2,b3,b4]
					coordinates.append(["TETRAHEDRAL","-","q is:",q])#= np.append(coordinates, [["TETRAHEDRAL","-","-","-"]], axis=0)
					atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('TETHDL',c+1,b1+1,b2+1,b3+1,b4+1,0,0)
					atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
					for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
						mcoordinates.append(fo[i,0:4])#=np.vstack((mcoordinates, fo[i,0:4]))
						coordinates.append(fo[i,0:4])#= np.vstack((coordinates, fo[i,0:4]))
						pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'TET','B',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
						pdbcoords.append(pdbtemp)#= np.vstack((pdbcoords, pdbtemp))
					coordinates.append([" "," "," "," "])#= np.append(coordinates, [[" "," "," "," "]], axis=0)
				elif   q <=0.85:
					q4fl.append(q)
					nlocal += 1
					n4fl += 1
					coordinates.append(["4-FOLD","-","q is:",q])#= np.append(coordinates, [["4-FOLD","-","-","-"]], axis=0)
					atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('4-FOLD',c+1,b1+1,b2+1,b3+1,b4+1,0,0)
					atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
					for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
						mcoordinates.append(fo[i,0:4])#=np.vstack((mcoordinates, fo[i,0:4]))
						coordinates.append(fo[i,0:4])#= np.vstack((coordinates, fo[i,0:4]))
						pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'4FL','C',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
						pdbcoords.append(pdbtemp)#= np.vstack((pdbcoords, pdbtemp))
					coordinates.append([" "," "," "," "])#= np.append(coordinates, [[" "," "," "," "]], axis=0)
				condition='true'
			if nbond==5:
				nlocal += 1
				n5fl += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2 + (cosb15 + 1/3)**2 + (cosb25 + 1/3)**2 + (cosb35 + 1/3)**2 + (cosb45 + 1/3)**2))
				#print(cosb12,cosb13,cosb14,cosb23,cosb24,cosb34,cosb15,cosb25,cosb35,cosb45)
				q5fl.append(q)
				aq.append(q)
				templist= [c,b1,b2,b3,b4,b5]
				#print("q n5:", q)
				coordinates.append(["5-FOLD","-","q is:",q])#= np.append(coordinates, [["5-FOLD","-","-","-"]], axis=0)
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('5-FOLD',c+1,b1+1,b2+1,b3+1,b4+1,b5+1,0)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4])#=np.vstack((mcoordinates, fo[i,0:4]))
					coordinates.append(fo[i,0:4])#= np.vstack((coordinates, fo[i,0:4]))
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'5FL','D',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)#= np.vstack((pdbcoords, pdbtemp))
				coordinates.append([" "," "," "," "])#= np.append(coordinates, [[" "," "," "," "]], axis=0)
				condition='true'
			if nbond==6:
				nlocal += 1
				noct += 1
				q=1-(0.375*((cosb12 + 1/3)**2 +(cosb13 + 1/3)**2 +(cosb14 + 1/3)**2 +(cosb23 + 1/3)**2 +(cosb24 + 1/3)**2 +(cosb34 + 1/3)**2 + (cosb15 + 1/3)**2 + (cosb25 + 1/3)**2 + (cosb35 + 1/3)**2 + (cosb45 + 1/3)**2 + (cosb16 + 1/3)**2 + (cosb26 + 1/3)**2 + (cosb36 + 1/3)**2 + (cosb46 + 1/3)**2 + (cosb56 + 1/3)**2))
				q6fl.append(q)
				aq.append(q)
				templist= [c,b1,b2,b3,b4,b5,b6]
				#print(templist)
				coordinates.append(["OCTAHEDRAL","-","q is:",q])#= np.append(coordinates, [["OCTAHEDRAL","-","-","-"]], axis=0)
				atomnumbertemp='{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('OCTHDL',c+1,b1+1,b2+1,b3+1,b4+1,b5+1,b6+1)
				atomnumbers.append(atomnumbertemp) # atom numbers, we will use it for pdos summation
				for  i in templist: #writes the coordinates in mcoord, coord and pdb-coord
					mcoordinates.append(fo[i,0:4])#=np.vstack((mcoordinates, fo[i,0:4]))
					coordinates.append(fo[i,0:4])#= np.vstack((coordinates, fo[i,0:4]))
					pdbtemp='{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',i+1,fo[i,0],'OCT','E',nlocal,fo[i,1],fo[i,2],fo[i,3],1.00,0.00,fo[i,0],"00")
					pdbcoords.append(pdbtemp)#= np.vstack((pdbcoords, pdbtemp))
				coordinates.append([" "," "," "," "])#= np.append(coordinates, [[" "," "," "," "]], axis=0)
				#print("q n6:", q)
				condition='true'
		else:
			continue
	Endf= "End of frame %d " % (N/Natom)
	mcoordinates.append([Endf," \r\n"," "," "])#=np.append(mcoordinates, [[Endf," \r\n"," "," "]], axis=0) #machine-readable coords, xyz coords of local orders
	stat= " \r\n Number of 3-FOLD, 4-FOLD, TETRAHEDRAL, 5-FOLD, OCTAHEDRAL: %d, %d, %d, %d, %d" % (n3fl, n4fl, ntet, n5fl, noct)
	locq= "\r\n Number of local orders are: %d " % (len(aq))
	coordinates.append([Endf,stat,locq,"\r\n"])#= np.append(coordinates, [[Endf,stat,locq,"\r\n"]], axis=0) # human-readable coords .
	localstats.append([n3fl,n4fl,ntet,n5fl,noct,nlocal])#= np.append(localstats, [[n3fl,n4fl,ntet,n5fl,noct,nlocal]],axis=0)
	pdbcoords.append([Endf])#= np.append(pdbcoords, [[Endf]], axis=0)
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
	print("End of frame %d " % (N/Natom))

#print (finalq, '\n' ,len(finalq[:]))
#afinalq= np.matrix(qframe)
#Tafinalq= qframe.T
#print(Tafinalq)
#qframe=[round(x,2) for x in qframe]
np.savetxt('output-q3fl.txt', q3flframe, delimiter=',', fmt="%s")
np.savetxt('output-q4fl.txt', q4flframe, delimiter=',', fmt="%s")
np.savetxt('output-qtet.txt', qtetframe, delimiter=',', fmt="%s")
np.savetxt('output-q5fl.txt', q5flframe, delimiter=',', fmt="%s")
np.savetxt('output-q6fl.txt', q6flframe, delimiter=',', fmt="%s")
np.savetxt('output-q-total.txt', qframe, delimiter=',', fmt="%s")
np.savetxt('output-human-readable-coords.txt', coordinates, delimiter=' ', fmt="%s")
np.savetxt('output-atom-number.txt',atomnumbers,delimiter=' ',fmt='%s')
np.savetxt('output-xyz-coords.txt', mcoordinates, delimiter=' ', fmt="%s")
np.savetxt('output-pdb-coords-m.txt', pdbcoords, fmt="%s")
np.savetxt('output-local-statistics.txt', localstats, fmt="%s")
np.savetxt('output-angle-distribution.txt', tetatot, fmt="%s")
print("Done!")
#np.savetxt('out field O- TiOx 0.1.txt',sorted_list,newline='\r\n') # baraye khate jadid bayad benevisi \r\n
