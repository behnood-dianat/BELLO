import numpy as np

file= np.loadtxt('output-human-readable-coords.txt', dtype='str', usecols = (0,1))

lenght= len(file)
nframe=0
nGe= 0
nTe= 0

n3FOLD=0
n4FOLD=0
nTET=0
n5FOLD=0
nOCT=0

n3FOLDGe=0
n4FOLDGe=0
nTETGe=0
n5FOLDGe=0
nOCTGe=0

n3FOLDTe=0
n4FOLDTe=0
nTETTe=0
n5FOLDTe=0
nOCTTe=0

Ge_3FL_3G0T=0
Ge_3FL_2G1T=0
Ge_3FL_1G2T=0
Ge_3FL_0G3T=0

Ge_4FL_4G0T=0
Ge_4FL_3G1T=0
Ge_4FL_2G2T=0
Ge_4FL_1G3T=0
Ge_4FL_0G4T=0

Ge_TET_4G0T=0
Ge_TET_3G1T=0
Ge_TET_2G2T=0
Ge_TET_1G3T=0
Ge_TET_0G4T=0

Ge_5FL_5G0T=0
Ge_5FL_4G1T=0
Ge_5FL_3G2T=0
Ge_5FL_2G3T=0
Ge_5FL_1G4T=0
Ge_5FL_0G5T=0

Ge_OCT_6G0T=0
Ge_OCT_5G1T=0
Ge_OCT_4G2T=0
Ge_OCT_3G3T=0
Ge_OCT_2G4T=0
Ge_OCT_1G5T=0
Ge_OCT_0G6T=0

Te_3FL_3G0T=0
Te_3FL_2G1T=0
Te_3FL_1G2T=0
Te_3FL_0G3T=0

Te_4FL_4G0T=0
Te_4FL_3G1T=0
Te_4FL_2G2T=0
Te_4FL_1G3T=0
Te_4FL_0G4T=0

Te_TET_4G0T=0
Te_TET_3G1T=0
Te_TET_2G2T=0
Te_TET_1G3T=0
Te_TET_0G4T=0

Te_5FL_5G0T=0
Te_5FL_4G1T=0
Te_5FL_3G2T=0
Te_5FL_2G3T=0
Te_5FL_1G4T=0
Te_5FL_0G5T=0

Te_OCT_6G0T=0
Te_OCT_5G1T=0
Te_OCT_4G2T=0
Te_OCT_3G3T=0
Te_OCT_2G4T=0
Te_OCT_1G5T=0
Te_OCT_0G6T=0

results= np.empty((0,66), int)
results=np.append(results,[["Ge","pcnt3foldGe","pcntGe_3FL_3G0T","pcntGe_3FL_2G1T","pcntGe_3FL_1G2T","pcntGe_3FL_0G3T","pcnt4foldGe","pcntGe_4FL_4G0T","pcntGe_4FL_3G1T","pcntGe_4FL_2G2T","pcntGe_4FL_1G3T","pcntGe_4FL_0G4T","pcnttetGe","pcntGe_TET_4G0T","pcntGe_TET_3G1T","pcntGe_TET_2G2T","pcntGe_TET_1G3T","pcntGe_TET_0G4T","pcnt5foldGe","pcntGe_5FL_5G0T","pcntGe_5FL_4G1T","pcntGe_5FL_3G2T","pcntGe_5FL_2G3T","pcntGe_5FL_1G4T","pcntGe_5FL_0G5T","pcntoctGe","pcntGe_OCT_6G0T","pcntGe_OCT_5G1T","pcntGe_OCT_4G2T","pcntGe_OCT_3G3T","pcntGe_OCT_2G4T","pcntGe_OCT_1G5T","pcntGe_OCT_0G6T","Te","pcnt3foldTe","pcntTe_3FL_3G0T","pcntTe_3FL_2G1T","pcntTe_3FL_1G2T","pcntTe_3FL_0G3T","pcnt4foldTe","pcntTe_4FL_4G0T","pcntTe_4FL_3G1T","pcntTe_4FL_2G2T","pcntTe_4FL_1G3T","pcntTe_4FL_0G4T","pcnttetTe","pcntTe_TET_4G0T","pcntTe_TET_3G1T","pcntTe_TET_2G2T","pcntTe_TET_1G3T","pcntTe_TET_0G4T","pcnt5foldTe","pcntTe_5FL_5G0T","pcntTe_5FL_4G1T","pcntTe_5FL_3G2T","pcntTe_5FL_2G3T","pcntTe_5FL_1G4T","pcntTe_5FL_0G5T","pcntoctTe","pcntTe_OCT_6G0T","pcntTe_OCT_5G1T","pcntTe_OCT_4G2T","pcntTe_OCT_3G3T","pcntTe_OCT_2G4T","pcntTe_OCT_1G5T","pcntTe_OCT_0G6T"]],axis=0)

for i in range(0,lenght):
    nGe = 0
    nTe = 0
    zerofilter=0
    if file[i,0] == 'End':
        #Ge percentages
        if (n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe) !=0:
            pcnt3foldGe= (n3FOLDGe*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcnt4foldGe= (n4FOLDGe*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcnttetGe= (nTETGe*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcnt5foldGe= (n5FOLDGe*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntoctGe= (nOCTGe*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)

            pcntGe_3FL_3G0T= (Ge_3FL_3G0T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_3FL_2G1T= (Ge_3FL_2G1T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_3FL_1G2T= (Ge_3FL_1G2T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_3FL_0G3T= (Ge_3FL_0G3T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)

            pcntGe_4FL_4G0T= (Ge_4FL_4G0T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_4FL_3G1T= (Ge_4FL_3G1T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_4FL_2G2T= (Ge_4FL_2G2T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_4FL_1G3T= (Ge_4FL_1G3T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_4FL_0G4T= (Ge_4FL_0G4T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)

            pcntGe_TET_4G0T= (Ge_TET_4G0T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_TET_3G1T= (Ge_TET_3G1T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_TET_2G2T= (Ge_TET_2G2T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_TET_1G3T= (Ge_TET_1G3T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_TET_0G4T= (Ge_TET_0G4T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)

            pcntGe_5FL_5G0T= (Ge_5FL_5G0T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_5FL_4G1T= (Ge_5FL_4G1T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_5FL_3G2T= (Ge_5FL_3G2T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_5FL_2G3T= (Ge_5FL_2G3T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_5FL_1G4T= (Ge_5FL_1G4T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_5FL_0G5T= (Ge_5FL_0G5T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)

            pcntGe_OCT_6G0T= (Ge_OCT_6G0T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_5G1T= (Ge_OCT_5G1T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_4G2T= (Ge_OCT_4G2T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_3G3T= (Ge_OCT_3G3T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_2G4T= (Ge_OCT_2G4T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_1G5T= (Ge_OCT_1G5T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
            pcntGe_OCT_0G6T= (Ge_OCT_0G6T*100)/(n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe)
        elif (n3FOLDGe+n4FOLDGe+nTETGe+n5FOLDGe+nOCTGe) ==0:
            pcnt3foldGe=pcntGe_3FL_3G0T=pcntGe_3FL_2G1T=pcntGe_3FL_1G2T=pcntGe_3FL_0G3T=pcnt4foldGe=pcntGe_4FL_4G0T=pcntGe_4FL_3G1T=pcntGe_4FL_2G2T=pcntGe_4FL_1G3T=pcntGe_4FL_0G4T=pcnttetGe=pcntGe_TET_4G0T=pcntGe_TET_3G1T=pcntGe_TET_2G2T=pcntGe_TET_1G3T=pcntGe_TET_0G4T=pcnt5foldGe=pcntGe_5FL_5G0T=pcntGe_5FL_4G1T=pcntGe_5FL_3G2T=pcntGe_5FL_2G3T=pcntGe_5FL_1G4T=pcntGe_5FL_0G5T=pcntoctGe=pcntGe_OCT_6G0T=pcntGe_OCT_5G1T=pcntGe_OCT_4G2T=pcntGe_OCT_3G3T=pcntGe_OCT_2G4T=pcntGe_OCT_1G5T=pcntGe_OCT_0G6T=0
            zerofilter += 1

        #Te percentages
        if (n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe) !=0:
            pcnt3foldTe= (n3FOLDTe*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcnt4foldTe= (n4FOLDTe*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcnttetTe= (nTETTe*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcnt5foldTe= (n5FOLDTe*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntoctTe= (nOCTTe*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)

            pcntTe_3FL_3G0T= (Te_3FL_3G0T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_3FL_2G1T= (Te_3FL_2G1T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_3FL_1G2T= (Te_3FL_1G2T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_3FL_0G3T= (Te_3FL_0G3T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)

            pcntTe_4FL_4G0T= (Te_4FL_4G0T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_4FL_3G1T= (Te_4FL_3G1T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_4FL_2G2T= (Te_4FL_2G2T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_4FL_1G3T= (Te_4FL_1G3T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_4FL_0G4T= (Te_4FL_0G4T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)

            pcntTe_TET_4G0T= (Te_TET_4G0T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_TET_3G1T= (Te_TET_3G1T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_TET_2G2T= (Te_TET_2G2T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_TET_1G3T= (Te_TET_1G3T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_TET_0G4T= (Te_TET_0G4T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)

            pcntTe_5FL_5G0T= (Te_5FL_5G0T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_5FL_4G1T= (Te_5FL_4G1T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_5FL_3G2T= (Te_5FL_3G2T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_5FL_2G3T= (Te_5FL_2G3T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_5FL_1G4T= (Te_5FL_1G4T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_5FL_0G5T= (Te_5FL_0G5T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)

            pcntTe_OCT_6G0T= (Te_OCT_6G0T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_5G1T= (Te_OCT_5G1T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_4G2T= (Te_OCT_4G2T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_3G3T= (Te_OCT_3G3T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_2G4T= (Te_OCT_2G4T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_1G5T= (Te_OCT_1G5T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
            pcntTe_OCT_0G6T= (Te_OCT_0G6T*100)/(n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe)
        elif (n3FOLDTe+n4FOLDTe+nTETTe+n5FOLDTe+nOCTTe) ==0:
            pcnt3foldTe=pcntTe_3FL_3G0T=pcntTe_3FL_2G1T=pcntTe_3FL_1G2T=pcntTe_3FL_0G3T=pcnt4foldTe=pcntTe_4FL_4G0T=pcntTe_4FL_3G1T=pcntTe_4FL_2G2T=pcntTe_4FL_1G3T=pcntTe_4FL_0G4T=pcnttetTe=pcntTe_TET_4G0T=pcntTe_TET_3G1T=pcntTe_TET_2G2T=pcntTe_TET_1G3T=pcntTe_TET_0G4T=pcnt5foldTe=pcntTe_5FL_5G0T=pcntTe_5FL_4G1T=pcntTe_5FL_3G2T=pcntTe_5FL_2G3T=pcntTe_5FL_1G4T=pcntTe_5FL_0G5T=pcntoctTe=pcntTe_OCT_6G0T=pcntTe_OCT_5G1T=pcntTe_OCT_4G2T=pcntTe_OCT_3G3T=pcntTe_OCT_2G4T=pcntTe_OCT_1G5T=pcntTe_OCT_0G6T=0
            zerofilter += 1

        if zerofilter == 2:
            continue
        results=np.append(results,[["Ge",pcnt3foldGe,pcntGe_3FL_3G0T,pcntGe_3FL_2G1T,pcntGe_3FL_1G2T,pcntGe_3FL_0G3T,pcnt4foldGe,pcntGe_4FL_4G0T,pcntGe_4FL_3G1T,pcntGe_4FL_2G2T,pcntGe_4FL_1G3T,pcntGe_4FL_0G4T,pcnttetGe,pcntGe_TET_4G0T,pcntGe_TET_3G1T,pcntGe_TET_2G2T,pcntGe_TET_1G3T,pcntGe_TET_0G4T,pcnt5foldGe,pcntGe_5FL_5G0T,pcntGe_5FL_4G1T,pcntGe_5FL_3G2T,pcntGe_5FL_2G3T,pcntGe_5FL_1G4T,pcntGe_5FL_0G5T,pcntoctGe,pcntGe_OCT_6G0T,pcntGe_OCT_5G1T,pcntGe_OCT_4G2T,pcntGe_OCT_3G3T,pcntGe_OCT_2G4T,pcntGe_OCT_1G5T,pcntGe_OCT_0G6T,"Te",pcnt3foldTe,pcntTe_3FL_3G0T,pcntTe_3FL_2G1T,pcntTe_3FL_1G2T,pcntTe_3FL_0G3T,pcnt4foldTe,pcntTe_4FL_4G0T,pcntTe_4FL_3G1T,pcntTe_4FL_2G2T,pcntTe_4FL_1G3T,pcntTe_4FL_0G4T,pcnttetTe,pcntTe_TET_4G0T,pcntTe_TET_3G1T,pcntTe_TET_2G2T,pcntTe_TET_1G3T,pcntTe_TET_0G4T,pcnt5foldTe,pcntTe_5FL_5G0T,pcntTe_5FL_4G1T,pcntTe_5FL_3G2T,pcntTe_5FL_2G3T,pcntTe_5FL_1G4T,pcntTe_5FL_0G5T,pcntoctTe,pcntTe_OCT_6G0T,pcntTe_OCT_5G1T,pcntTe_OCT_4G2T,pcntTe_OCT_3G3T,pcntTe_OCT_2G4T,pcntTe_OCT_1G5T,pcntTe_OCT_0G6T]],axis=0)

        nframe += 1

        n3FOLD=0
        n4FOLD=0
        nTET=0
        n5FOLD=0
        nOCT=0

        n3FOLDGe=0
        n4FOLDGe=0
        nTETGe=0
        n5FOLDGe=0
        nOCTGe=0

        n3FOLDTe=0
        n4FOLDTe=0
        nTETTe=0
        n5FOLDTe=0
        nOCTTe=0

        Ge_3FL_3G0T=0
        Ge_3FL_2G1T=0
        Ge_3FL_1G2T=0
        Ge_3FL_0G3T=0

        Ge_4FL_4G0T=0
        Ge_4FL_3G1T=0
        Ge_4FL_2G2T=0
        Ge_4FL_1G3T=0
        Ge_4FL_0G4T=0

        Ge_TET_4G0T=0
        Ge_TET_3G1T=0
        Ge_TET_2G2T=0
        Ge_TET_1G3T=0
        Ge_TET_0G4T=0

        Ge_5FL_5G0T=0
        Ge_5FL_4G1T=0
        Ge_5FL_3G2T=0
        Ge_5FL_2G3T=0
        Ge_5FL_1G4T=0
        Ge_5FL_0G5T=0

        Ge_OCT_6G0T=0
        Ge_OCT_5G1T=0
        Ge_OCT_4G2T=0
        Ge_OCT_3G3T=0
        Ge_OCT_2G4T=0
        Ge_OCT_1G5T=0
        Ge_OCT_0G6T=0

        Te_3FL_3G0T=0
        Te_3FL_2G1T=0
        Te_3FL_1G2T=0
        Te_3FL_0G3T=0

        Te_4FL_4G0T=0
        Te_4FL_3G1T=0
        Te_4FL_2G2T=0
        Te_4FL_1G3T=0
        Te_4FL_0G4T=0

        Te_TET_4G0T=0
        Te_TET_3G1T=0
        Te_TET_2G2T=0
        Te_TET_1G3T=0
        Te_TET_0G4T=0

        Te_5FL_5G0T=0
        Te_5FL_4G1T=0
        Te_5FL_3G2T=0
        Te_5FL_2G3T=0
        Te_5FL_1G4T=0
        Te_5FL_0G5T=0

        Te_OCT_6G0T=0
        Te_OCT_5G1T=0
        Te_OCT_4G2T=0
        Te_OCT_3G3T=0
        Te_OCT_2G4T=0
        Te_OCT_1G5T=0
        Te_OCT_0G6T=0
        #continue
    else:
        if file[i,0] =='3-FOLD':
            n3FOLD +=1
            if file[i+1,0] == 'Ge':
                n3FOLDGe +=1
                for j in range(i+2,i+5):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    elif file[j,0] == 'Te':
                        nTe += 1

                if nGe==3 and nTe==0:
                    Ge_3FL_3G0T += 1
                if nGe==2 and nTe==1:
                    Ge_3FL_2G1T += 1
                if nGe==1 and nTe==2:
                    Ge_3FL_1G2T += 1
                if nGe==0 and nTe==3:
                    Ge_3FL_0G3T += 1

            if file[i+1,0] == 'Te':
                n3FOLDTe +=1
                for j in range(i+2,i+5):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==3 and nTe==0:
                    Te_3FL_3G0T += 1
                if nGe==2 and nTe==1:
                    Te_3FL_2G1T += 1
                if nGe==1 and nTe==2:
                    Te_3FL_1G2T += 1
                if nGe==0 and nTe==3:
                    Te_3FL_0G3T += 1
                #---------------------------------------------------------------
        if file[i,0] =='4-FOLD':
            n4FOLD +=1
            if file[i+1,0] == 'Ge':
                n4FOLDGe +=1
                for j in range(i+2,i+6):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==4 and nTe==0:
                    Ge_4FL_4G0T += 1
                if nGe==3 and nTe==1:
                    Ge_4FL_3G1T += 1
                if nGe==2 and nTe==2:
                    Ge_4FL_2G2T += 1
                if nGe==1 and nTe==3:
                    Ge_4FL_1G3T += 1
                if nGe==0 and nTe==4:
                    Ge_4FL_0G4T += 1

            if file[i+1,0] == 'Te':
                n4FOLDTe +=1
                for j in range(i+2,i+6):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==4 and nTe==0:
                    Te_4FL_4G0T += 1
                if nGe==3 and nTe==1:
                    Te_4FL_3G1T += 1
                if nGe==2 and nTe==2:
                    Te_4FL_2G2T += 1
                if nGe==1 and nTe==3:
                    Te_4FL_1G3T += 1
                if nGe==0 and nTe==4:
                    Te_4FL_0G4T += 1

        if file[i,0] =='TETRAHEDRAL':
            nTET +=1
            if file[i+1,0] == 'Ge':
                nTETGe +=1
                for j in range(i+2,i+6):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==4 and nTe==0:
                    Ge_TET_4G0T += 1
                if nGe==3 and nTe==1:
                    Ge_TET_3G1T += 1
                if nGe==2 and nTe==2:
                    Ge_TET_2G2T += 1
                if nGe==1 and nTe==3:
                    Ge_TET_1G3T += 1
                if nGe==0 and nTe==4:
                    Ge_TET_0G4T += 1

            if file[i+1,0] == 'Te':
                nTETTe +=1
                for j in range(i+2,i+6):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==4 and nTe==0:
                    Te_TET_4G0T += 1
                if nGe==3 and nTe==1:
                    Te_TET_3G1T += 1
                if nGe==2 and nTe==2:
                    Te_TET_2G2T += 1
                if nGe==1 and nTe==3:
                    Te_TET_1G3T += 1
                if nGe==0 and nTe==4:
                    Te_TET_0G4T += 1

        if file[i,0] =='5-FOLD':
            n5FOLD +=1
            if file[i+1,0] == 'Ge':
                n5FOLDGe +=1
                for j in range(i+2,i+7):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==5 and nTe==0:
                    Ge_5FL_5G0T += 1
                if nGe==4 and nTe==1:
                    Ge_5FL_4G1T += 1
                if nGe==3 and nTe==2:
                    Ge_5FL_3G2T += 1
                if nGe==2 and nTe==3:
                    Ge_5FL_2G3T += 1
                if nGe==1 and nTe==4:
                    Ge_5FL_1G4T += 1
                if nGe==0 and nTe==5:
                    Ge_5FL_0G5T += 1

            if file[i+1,0] == 'Te':
                n5FOLDTe +=1
                for j in range(i+2,i+7):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==5 and nTe==0:
                    Te_5FL_5G0T += 1
                if nGe==4 and nTe==1:
                    Te_5FL_4G1T += 1
                if nGe==3 and nTe==2:
                    Te_5FL_3G2T += 1
                if nGe==2 and nTe==3:
                    Te_5FL_2G3T += 1
                if nGe==1 and nTe==4:
                    Te_5FL_1G4T += 1
                if nGe==0 and nTe==5:
                    Te_5FL_0G5T += 1

        if file[i,0] =='OCTAHEDRAL':
            nOCT +=1
            if file[i+1,0] == 'Ge':
                nOCTGe +=1
                for j in range(i+2,i+8):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==6 and nTe==0:
                    Ge_OCT_6G0T += 1
                if nGe==5 and nTe==1:
                    Ge_OCT_5G1T += 1
                if nGe==4 and nTe==2:
                    Ge_OCT_4G2T += 1
                if nGe==3 and nTe==3:
                    Ge_OCT_3G3T += 1
                if nGe==2 and nTe==4:
                    Ge_OCT_2G4T += 1
                if nGe==1 and nTe==5:
                    Ge_OCT_1G5T += 1
                if nGe==0 and nTe==6:
                    Ge_OCT_0G6T += 1

            if file[i+1,0] == 'Te':
                nOCTTe +=1
                for j in range(i+2,i+8):
                    if file[j,0] == 'Ge':
                        nGe += 1
                    if file[j,0] == 'Te':
                        nTe += 1

                if nGe==6 and nTe==0:
                    Te_OCT_6G0T += 1
                if nGe==5 and nTe==1:
                    Te_OCT_5G1T += 1
                if nGe==4 and nTe==2:
                    Te_OCT_4G2T += 1
                if nGe==3 and nTe==3:
                    Te_OCT_3G3T += 1
                if nGe==2 and nTe==4:
                    Te_OCT_2G4T += 1
                if nGe==1 and nTe==5:
                    Te_OCT_1G5T += 1
                if nGe==0 and nTe==6:
                    Te_OCT_0G6T += 1

        else:
            continue

np.savetxt('COORDINATION-STATS.txt', results, delimiter=' ', fmt="%s")
