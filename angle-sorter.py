import numpy as np

file= np.loadtxt('output-angle-distribution.txt', dtype='str', usecols = (0,1))

lenght=len(file)
tetete=[]
tetege=[]
tegete=[]
getete=[]
gegege=[]
getege=[]
gegete=[]
tegege=[]
print("lenght is:",lenght)
for i in range(0,lenght):
    if file[i,1]=='Te-Te-Te':
        tetete.append(file[i])
    if file[i,1]=='Te-Te-Ge':
        tetege.append(file[i])
    if file[i,1]=='Te-Ge-Te':
        tegete.append(file[i])
    if file[i,1]=='Ge-Te-Te':
        getete.append(file[i])
    if file[i,1]=='Ge-Ge-Ge':
        gegege.append(file[i])
    if file[i,1]=='Ge-Te-Ge':
        getege.append(file[i])
    if file[i,1]=='Te-Ge-Ge':
        tegege.append(file[i])
    if file[i,1]=='Ge-Ge-Te':
        gegete.append(file[i])

np.savetxt('tetete.txt', tetete, fmt="%s")
np.savetxt('tetege.txt', tetege, fmt="%s")
np.savetxt('tegete.txt', tegete, fmt="%s")
np.savetxt('getete.txt', getete, fmt="%s")
np.savetxt('gegege.txt', gegege, fmt="%s")
np.savetxt('getege.txt', getege, fmt="%s")
np.savetxt('gegete.txt', gegete, fmt="%s")
np.savetxt('tegege.txt', tegege, fmt="%s")
