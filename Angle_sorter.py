import shutil
import numpy as np
import os, pathlib, stat
import matplotlib.pyplot as plt
from itertools import product
import scipy.stats as stats

def folder_overwrite(sorted_folder):
    if os.path.exists(sorted_folder):
        os.chmod(sorted_folder, stat.S_IWRITE)
        shutil.rmtree(sorted_folder)
    os.mkdir(sorted_folder)

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

def sorter(elements):
    file= np.loadtxt('output-angle-distribution.txt',delimiter=' ', dtype="f,U8", usecols = (0,1))

    lenght=len(file)

    folder_path = pathlib.PurePosixPath(os.getcwd())
    sorted_folder = pathlib.Path(folder_path, 'Sorted_angles')
    folder_overwrite(sorted_folder)

    names=[]
    file_names=[]
    plotting_list=[]

    def func_element(elements):
        sorted_list = []
        for i in elements:
            for j in elements:
                for k in elements:
                    name=i+'-'+j+'-'+k
                    for x in range(0,lenght):
                        if file[x][1]==name:
                            sorted_list.append(list(file[x]))
                    save_name = name + '.txt'
                    names.append(name)
                    file_names.append(save_name)
                    np.savetxt(pathlib.PurePosixPath(sorted_folder, save_name), sorted_list, fmt="%s")
                    plotting_list.append(sorted_list)
                    sorted_list=[]


    func_element(elements)
    print('Angle sorting is done!')
    bins = np.linspace(0, 180, 45, dtype='i')
    for i in range(0,len(file_names)):
       plt.figure(i)
       data = np.loadtxt(pathlib.PurePosixPath(sorted_folder, file_names[i]), dtype="f", usecols=0)
       n, x, patches = plt.hist(x=data,bins=bins,label=names[i],
                                              edgecolor='black',density=True)
       cm = plt.cm.get_cmap('viridis')
       col = (n - n.min()) / (n.max() - n.min())
       for c, p in zip(col, patches):
           plt.setp(p, 'facecolor', cm(c))
       density = stats.gaussian_kde(data)
       plt.plot(x, density(x),label=('Fit %s' % names[i]))
       plt.legend()

 #  -----------------------------------------------------------------------------------------------
 #  ---- This is for sub-plotting, it's best suited for 2-element plots since more than 2 elements
 #  ---- will yield a lot of combinations----------------------------------------------------------
 #  -----------------------------------------------------------------------------------------------
 #  def subplot_plotting(elements):
 #      elements_combinations=len(list(product(elements, repeat=3)))
 #      if elements_combinations==8:
 #          tuples=[]
 #          shape=(2,4)
 #          for i in range(0, 2):
 #              for j in range(0, 4):
 #                  tuples.append((i, j))
 #      elif elements_combinations==27:
 #          tuples=[]
 #          shape=(3,9)
 #          for i in range(0, 3):
 #              for j in range(0, 9):
 #                  tuples.append((i, j))
 #      fig, ax = plt.subplots(shape[0],shape[1])
 #      for i in range(0, len(file_names)):
 #          data = np.loadtxt(pathlib.PurePosixPath(sorted_folder, file_names[i]), dtype="f", usecols=0)
 #          n, x, patches = ax[tuples[i]].hist(x=data, bins=bins, label=names[i],
 #                                   edgecolor='black', density=True)
 #          cm = plt.cm.get_cmap('viridis')
 #          col = (n - n.min()) / (n.max() - n.min())
 #          for c, p in zip(col, patches):
 #              plt.setp(p, 'facecolor', cm(c))
 #          density = stats.gaussian_kde(data)
 #          ax[tuples[i]].plot(x, density(x), label=('Fit %s' % names[i]),color='black')
 #          ax[tuples[i]].legend()
 #          subplot_plotting(elements)

    plt.show()