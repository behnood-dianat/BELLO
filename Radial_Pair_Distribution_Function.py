import numpy as np
import math as mt
import pandas as pd
import matplotlib.pyplot as plt
import PySimpleGUI as sg


def RDF(url,celldmx_raw,celldmy_raw,celldmz_raw,dr,rmax,first_element,second_element):
	print("|-----------------------------------------------------|\r\n"
		  "|----------------------B.E.L.L.O----------------------|\r\n"
		  "|---------Bond Element Lattice Locality Order---------|\r\n"
		  "|-----------------------------------------------------|\r\n"
		  "|----------Radial Pair Distribution Function----------|")
	f = pd.read_fwf(url,header=None)
	f= f.fillna("x")
	file = np.array(f)
	Natom=int(f[0][0])

	#------------------------------------
	#--------------Inputs----------------
	#------------------------------------
	celldmx=float(celldmx_raw)
	celldmy=float(celldmy_raw)
	celldmz=float(celldmz_raw)
	dr=float(dr)
	rmax=float(rmax)

	fat=first_element
	sat=second_element

	lfile=len(file)
	print("File lenght is: ",lfile)
	#------------------------------------
	#-----variables and empy lists-------
	#------------------------------------
	fa= np.empty((0,4), int)
	zerocenter=[0,0,0]
	templist=[]
	div=10
	intervals=1/dr
	Nframes= round(len(f)/(Natom+2)) #number of frames
	print("Number of frames: ",Nframes)
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
	print("Boundry Condition is done! \r\nCalculating Radial Pair Distribution Function:")

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
	totalRDF=[]
	#-------------------------------------
	#---looping over the atoms of a frame-
	#-------------------------------------
	for N in range(0,len(fa),Natom):
	#-variables for each frame-----
		sg.one_line_progress_meter('Frames', current_value=N, max_value=len(fa), orientation='h')
		rlist=[]
		fo=[]
		fi=[]
		natm1=0
		natm2=0
		fo=np.copy(fa[N:N+Natom])
		lfo=len(fo)
		g=[]
		d=10
		for l in range(lfo):
			if fo[l,0]==fat:
				natm1 += 1
			if fo[l,0]==sat:
				natm2 += 1


		for c in range (0,lfo):
			if fo[c,0] == fat:
				for b1 in range(0,lfo):
					if (c!=b1) & (fo[b1,0] == sat):
						a=fo[c,1:4]
						b=fo[b1,1:4]
						db1=distance(a,b)
						templist.append(db1)
						if  (db1 <= rmax) and (db1 <= celldmx/2) and (db1 <= celldmy/2) and (db1 <= celldmz/2): #maximum range condition
							rlist.append(db1)

		rmax=(max(rlist))
		nvals=len(rlist)
		n= int(mt.ceil(rmax/dr))
		rdf_hist = np.zeros(n)
		rr = np.arange(0+dr,rmax+dr,dr)
		r = np.around(rr,decimals=3)
		area = [4.*(mt.pi)*(r[i]**2) for i in range(n)]
		vol = [(4.0 / 3.0) * (mt.pi) * (r[i])**3 for i in range(n)]
		vol.insert(0,0)
		density=natm2/(celldmx*celldmy*celldmz)

		for i in range(n):
			for j in range(nvals):
				test = (round(rlist[j]*intervals)/intervals) # round the "distance" to the nearest interval(dr)
				if (test == r[i]):
					rdf_hist[i] += 1

		for i in range(n):
			#rdf_hist[i] /= area
			rdf_hist[i] /= natm1
			rdf_hist[i] /= density
			rdf_hist[i] /= vol[i]-vol[i-1]
		totalRDF.append(rdf_hist)

	rdf=sum(totalRDF)/Nframes
	print(len(fi))
	finalfile= np.column_stack((r,rdf))
#
	np.savetxt('RDF.txt', finalfile, delimiter=' ', fmt="%s")
	sg.one_line_progress_meter_cancel()
	fig = plt.figure(facecolor="none",figsize=(6.75, 5))
	ax = fig.add_subplot(111)
	fontsize = 14
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(fontsize)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(fontsize)
	plt.plot(r,rdf,label='radial distribution function',lw=2)
	plt.xlabel('radius (Angstroms)',fontsize=fontsize)
	plt.ylabel('radial pair distribution function g(r)',fontsize=fontsize)
	plt.show()