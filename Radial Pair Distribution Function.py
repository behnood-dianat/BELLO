import numpy as np
import math as mt
import pandas as pd
import matplotlib.pyplot as plt


print("|-----------------------------------------------------| \r\n|----------------------B.E.L.L.O----------------------| \r\n|---------Bond Element Lattice Locality Order---------|\r\n|-----------------------------------------------------|\r\n|----------Radial Pair Distribution Function----------|")

f = pd.read_fwf('aaa.xyz',header=None)
f= f.fillna("x") 
file = np.array(f)
Natom=int(f[0][0])
#emptylist=[" "," "," "," "]
celldm=20.644848284633
#celldm=float(input('Enter lattice constant (in Ang): '))
#fat=str(input('Name of the first atom: '))
#sat=str(input('Name of the second atom: '))
#limit=float(input('Maximum range: '))
fat="Ge"
sat="Te"
rmax=10
dr=0.1

div=10
intervals=1/dr


lfile=len(file)
print("File lenght is: ",lfile)

#------------------------------------
#-----variables and empy lists-------
#------------------------------------
fa= np.empty((0,4), int)
zerocenter=[0,0,0]
templist=[]

#------------------------------------
#------boundry condition-------------
#------------------------------------
for i in range(lfile):
	if (file[i,1] != 'x'):
		if file[i,1]>celldm or file[i,1]<0:
			file[i,1]= file[i,1]%celldm
		if file[i,2]>celldm or file[i,2]<0:
			file[i,2]= file[i,2]%celldm
		if file[i,3]>celldm or file[i,3]<0:
			file[i,3]= file[i,3]%celldm
		fa= np.vstack((fa, file[i]))
print("Boundry Condition is done! \r\nCalculating local orders:")


#-------------------------------------
#---looping over the atoms of a frame-
#-------------------------------------
for N in range(0,len(fa),Natom):
#-variables for each frame-----
	rlist=[]
	fo=[]
	fi=[]
	natm1=0
	natm2=0
	fo=np.copy(fa[N:N+Natom])
	#fi=[]
#	fi=np.copy(fa[N:N+Natom])
	lfo=len(fo)
	g=[]
	d=10 #d >= 2
	for l in range(lfo):
		if fo[l,0]==fat:
			natm1 += 1
		if fo[l,0]==sat:
			natm2 += 1
	def distance(a, b):
		dx = abs(a[0] - b[0])
		x = min(dx, abs(A - dx))
		
		dy = abs(a[1] - b[1])
		y = min(dy, abs(B - dy))
		
		dz = abs(a[2] - b[2])
		z = min(dz, abs(C - dz))

		return mt.sqrt(x**2 + y**2 + z**2)

	A=B=C=celldm

	for c in range (0,lfo):
		if fo[c,0] == fat:
			for b1 in range(0,lfo):
				if (c!=b1) & (fo[b1,0] == sat):
					a=fo[c,1:4]
					b=fo[b1,1:4]
					db1=distance(a,b)
					templist.append(db1)
					if  (db1 <= rmax) and (db1 <= celldm/2): #maximum range condition
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
	density=natm2/(celldm**3)

	for i in range(n):
		for j in range(nvals):
			test = (round(rlist[j]*intervals)/intervals) # round the "distance" to the nearest interval(dr)
			if (test == r[i]):
				rdf_hist[i] += 1 
				
	for i in range(n):
		#rdf_hist[i] /= area[i]
		rdf_hist[i] /= natm1
		rdf_hist[i] /= density 
		rdf_hist[i] /= vol[i]-vol[i-1]
		
print(len(fi))
finalfile= np.column_stack((r,rdf_hist))

np.savetxt('RPDF.txt', finalfile, delimiter=',', fmt="%s")
np.savetxt('RPDF-xyz.txt', fi, delimiter=' ', fmt="%s")

plt.plot(r,rdf_hist,label='radial distribution function')
plt.xlabel('radius (Angstroms)')
plt.ylabel('radial pair distribution function g(r)')
plt.show()