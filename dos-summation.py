import numpy as np
from numpy import sum
import os
import glob
import string

atn= np.loadtxt('output-atom-number.txt', delimiter=" ",dtype={'names': ('atom', 'center', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6'), 'formats':('U6', 'float', 'float', 'float', 'float', 'float', 'float', 'float')})


filename= "Ge50Se50_100atoms.pdos_atm#" + "1"+ "(" + "*"         #|    --------------------------------------------
file=glob.glob(filename)                                   #|    Since all the energies are of the same range,  
pdos= np.genfromtxt(file[0], delimiter="  ",dtype='f')     #|    we gather it before going into the loops
energy=pdos[:,0]                                           #|    which saves a lot of computational power
a3folds=np.zeros([len(pdos),5])          # all of 3Folds in one file per (E,S,P,D,T(total=s+p+d))
a3folds[:,0] = energy[:]                 # put the energy in first coloumn
a4folds=np.zeros([len(pdos),5])
a4folds[:,0] = energy[:]
atet=np.zeros([len(pdos),5])
atet[:,0] = energy[:]
a5folds=np.zeros([len(pdos),5])
a5folds[:,0] = energy[:]
aoct=np.zeros([len(pdos),5])
aoct[:,0] = energy[:]
a3456=np.zeros([len(pdos),6])            # energy per total (s+p+d) of 3fold/4fold/tetrahedral/5fold/octahedral
a3456[:,0] = energy[:]

title="E"+" "+"S"+" "+"P"+" "+"D"+" "+"T"
title2="3FOLD"+" "+"4FOLD"+" "+"TET"+" "+"5FOLD"+" "+"OCT"

for i in range (len(atn)): # each localoty (3fold, 4fold,...)
    print(i)
    sc=[]     # coloumn for s orbital
    pc=[]     # coloumn for p orbital
    dc=[]     # coloumn for d orbital
    tspd=[]    # = sc + pc + dc
    for j in range(1,8): # each atom of the locality
        if atn[i][j] != 0: # 0 means no atoms
            filename= "Ge50Se50_100atoms.pdos_atm#"+str(int(atn[i][j]))+"("+"*" 
            file=glob.glob(filename)
            if len(file) == 0:
                    print("Missing file")
            for k in file: # each atom has multiple orbitals (S, P,...)
                pdos= np.genfromtxt(k, delimiter="  ",dtype='f')
                if len(pdos[0]) == 3:
                    sc.append(pdos[:,2])
                    tempsc=np.array(sc)
                    sc.clear()
                    switch= [0]
                elif len(pdos[0]) == 5:
                    pc.append(pdos[:,2:5])
                    temppc=np.array(pc)
                    pc.clear()
                    switch= [1]
                elif len(pdos[0]) == 7:
                    dc.append(pdos[:,2:7])
                    tempdc=np.array(dc)
                    dc.clear()
                    switch= [2]
                else:
                    pdos=np.empty(pdos.shape)
                    print("wrong coloumns")
                if atn[i][0] == '3-FOLD':
                    if switch == [0]:
                        a3folds[:,1] = a3folds[:,1] + tempsc[0,:]
                        del(tempsc,switch)
                    elif switch == [1]:
                        a3folds[:,2] = a3folds[:,2] + temppc[0,:,0] + temppc[0,:,1] + temppc[0,:,2]
                        del(temppc,switch)
                    elif switch == [2]:
                        a3folds[:,3] = a3folds[:,3] + tempdc[0,:,0] + tempdc[0,:,1] + tempdc[0,:,2] + tempdc[0,:,3] + tempdc[0,:,4]
                        del(tempdc,switch)
                elif atn[i][0] == '4-FOLD':
                    if switch == [0]:
                        a4folds[:,1] = a4folds[:,1] + tempsc[0,:]
                        del(tempsc,switch)
                    elif switch == [1]:
                        a4folds[:,2] = a4folds[:,2] + temppc[0,:,0] + temppc[0,:,1] + temppc[0,:,2]
                        del(temppc,switch)
                    elif switch == [3]:
                        a4folds[:,3] = a4folds[:,3] + tempdc[0,:,0] + tempdc[0,:,1] + tempdc[0,:,2] + tempdc[0,:,3] + tempdc[0,:,4]
                        del(tempdc,switch)
                elif atn[i][0] == 'TETHDL':
                    if switch == [0]:
                        atet[:,1] = atet[:,1] + tempsc[0,:]
                        del(tempsc,switch)
                    elif switch == [1]:
                        atet[:,2] = atet[:,2] + temppc[0,:,0] + temppc[0,:,1] + temppc[0,:,2]
                        del(temppc,switch)
                    elif switch == [3]:
                        atet[:,3] = atet[:,3] + tempdc[0,:,0] + tempdc[0,:,1] + tempdc[0,:,2] + tempdc[0,:,3] + tempdc[0,:,4]
                        del(tempdc,switch)
                elif atn[i][0] == '5-FOLD':
                    if switch == [0]:
                        a5folds[:,1] = a5folds[:,1] + tempsc[0,:]
                        del(tempsc,switch)
                    elif switch == [1]:
                        a5folds[:,2] = a5folds[:,2] + temppc[0,:,0] + temppc[0,:,1] + temppc[0,:,2]
                        del(temppc,switch)
                    elif switch == [3]:
                        a5folds[:,3] = a5folds[:,3] + tempdc[0,:,0] + tempdc[0,:,1] + tempdc[0,:,2] + tempdc[0,:,3] + tempdc[0,:,4]
                        del(tempdc,switch)
                elif atn[i][0] == 'OCTHDL':
                    if switch == [0]:
                        aoct[:,1] = aoct[:,1] + tempsc[0,:]
                        del(tempsc,switch)
                    elif switch == [1]:
                        aoct[:,2] = aoct[:,2] + temppc[0,:,0] + temppc[0,:,1] + temppc[0,:,2]
                        del(temppc,switch)
                    elif switch == [3]:
                        aoct[:,3] = aoct[:,3] + tempdc[0,:,0] + tempdc[0,:,1] + tempdc[0,:,2] + tempdc[0,:,3] + tempdc[0,:,4]
                        del(tempdc,switch)
                del(pdos)
a3folds[:,4] = a3folds[:,1] + a3folds[:,2] + a3folds[:,3]
a4folds[:,4] = a4folds[:,1] + a4folds[:,2] + a4folds[:,3]                
a5folds[:,4] = a5folds[:,1] + a5folds[:,2] + a5folds[:,3]
atet[:,4] = atet[:,1] + atet[:,2] + atet[:,3]
aoct[:,4] = aoct[:,1] + aoct[:,2] + aoct[:,3]
a3456[:,1]=a3folds[:,4]
a3456[:,2]=a4folds[:,4]
a3456[:,3]=atet[:,4]
a3456[:,4]=a5folds[:,4]
a3456[:,5]=aoct[:,4]
np.savetxt('All-3Folds.txt', a3folds, delimiter=' ', fmt="%s", header=title)
np.savetxt('All-4Folds.txt', a4folds, delimiter=' ', fmt="%s", header=title)
np.savetxt('All-5Folds.txt', a5folds, delimiter=' ', fmt="%s", header=title)
np.savetxt('All-Tetrahedrals.txt', atet, delimiter=' ', fmt="%s", header=title)
np.savetxt('All-Octahedrals.txt', aoct, delimiter=' ', fmt="%s", header=title)
np.savetxt('All-Localities.txt', a3456, delimiter=' ', fmt="%s", header=title2)
print("Done!")