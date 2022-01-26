import numpy as np
import pandas as pd
from collections import Counter
from itertools import combinations_with_replacement as cr
import seaborn as sns
import matplotlib.pyplot as plt

def coordination_heatmap(*elements_raw):
	file= np.loadtxt('output-human-readable-coords.txt', dtype='str', usecols = (0,1))
	Nframe=int(file[0,1])

	def merger(input_list):
		if isinstance(input_list,list):
			input_list.sort() # for lists
		else:
			sorted(input_list)# for tuples
		output=''.join([f'{key}{value}' for key, value in Counter(input_list).items()])
		return(output)

	def string_attach(input_list):
		output=[]
		for i in input_list:
			output.append(merger(i))
		return(output)

	def column_lenght(elements,fold):
		output=len(list(cr(elements,fold)))
		return(output)

	# Set Dataframe indexes-------------------------------------------------
	def indexer(grid):
		fold3_index=string_attach(list(cr(elements,3)))
		grid.loc[:column_lenght(elements,3),'3_fold_index']=fold3_index
		fold4_index=string_attach(list(cr(elements,4)))
		grid.loc[:column_lenght(elements,4),'4_fold_index']	 =fold4_index
		grid.loc[:column_lenght(elements,4),'Tetrahedral_index']=fold4_index
		fold5_index=string_attach(list(cr(elements,5)))
		grid.loc[:column_lenght(elements,5),'5_fold_index']=fold5_index
		fold6_index=string_attach(list(cr(elements,6)))
		grid.loc[:column_lenght(elements,6),'6_fold_index']=fold6_index
		return(grid)
	#-----------------------------------------------------------------------
	elements=list(elements_raw)
	elements = [item for sublist in elements for item in sublist]
	elements.sort()

	columns=['3_fold_index','3_fold','4_fold_index','4_fold','Tetrahedral_index','Tetrahedral','5_fold_index','5_fold','6_fold_index','6_fold']
	lenght=column_lenght(elements,6)
	grid=pd.DataFrame(0,index=np.linspace(1,lenght,lenght,dtype='i'),columns=columns)
	grid=indexer(grid)



	#-Coordination Heatmap irrelevant of the center atom specie---------------------------------------------
	for i in range(0,len(file)):
		if file[i,0] =='3-FOLD':
			temp_element=[file[i+2,0],file[i+3,0],file[i+4,0]]
			for j in range(1,(column_lenght(elements,3))):
				if merger(temp_element) == grid.loc[j,'3_fold_index']:
					grid.loc[j,'3_fold'] += 1
		elif file[i,0] =='4-FOLD':
			temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0]]
			for j in range(1,(column_lenght(elements,4))):
				if merger(temp_element) == grid.loc[j,'4_fold_index']:
					grid.loc[j,'4_fold'] += 1
		elif file[i,0] =='TETRAHEDRAL':
			temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0]]
			for j in range(1,(column_lenght(elements,4))):
				if merger(temp_element) == grid.loc[j,'Tetrahedral_index']:
					grid.loc[j,'Tetrahedral'] += 1
		elif file[i,0] =='5-FOLD':
			temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0],file[i+6,0]]
			for j in range(1,(column_lenght(elements,5))):
				if merger(temp_element) == grid.loc[j,'5_fold_index']:
					grid.loc[j,'5_fold'] += 1
		elif file[i,0] =='OCTAHEDRAL':
			temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0],file[i+6,0],file[i+7,0]]
			for j in range(1,(column_lenght(elements,6))):
				if merger(temp_element) == grid.loc[j,'6_fold_index']:
					grid.loc[j,'6_fold'] += 1
	grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')]=(grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')])/Nframe
	plt.figure()
	sns.heatmap(grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')],
		annot=(grid.loc[:,('3_fold_index','4_fold_index','Tetrahedral_index','5_fold_index','6_fold_index')]),
		linewidths=.5,cmap="viridis",fmt='',square=False, yticklabels=False)
	plt.title("Coordination Heatmap irrelevant of the center atom specie")
	#--------------------------------------------------------------------------------------------------------


	#-Coordination Heatmap with elements as center atom-----------------------------------------------------
	for x in elements:
		grid=indexer(pd.DataFrame(0,index=np.linspace(1,lenght,lenght,dtype='i'),columns=columns))
		for i in range(0,len(file)):
			if file[i,0] =='3-FOLD' and file[i+1,0]==x:
				temp_element=[file[i+2,0],file[i+3,0],file[i+4,0]]
				for j in range(1,(column_lenght(elements,3))):
					if merger(temp_element) == grid.loc[j,'3_fold_index']:
						grid.loc[j,'3_fold'] += 1
			elif file[i,0] =='4-FOLD' and file[i+1,0]==x:
				temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0]]
				for j in range(1,(column_lenght(elements,4))):
					if merger(temp_element) == grid.loc[j,'4_fold_index']:
						grid.loc[j,'4_fold'] += 1
			elif file[i,0] =='TETRAHEDRAL' and file[i+1,0]==x:
				temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0]]
				for j in range(1,(column_lenght(elements,4))):
					if merger(temp_element) == grid.loc[j,'Tetrahedral_index']:
						grid.loc[j,'Tetrahedral'] += 1
			elif file[i,0] =='5-FOLD' and file[i+1,0]==x:
				temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0],file[i+6,0]]
				for j in range(1,(column_lenght(elements,5))):
					if merger(temp_element) == grid.loc[j,'5_fold_index']:
						grid.loc[j,'5_fold'] += 1
			elif file[i,0] =='OCTAHEDRAL' and file[i+1,0]==x:
				temp_element=[file[i+2,0],file[i+3,0],file[i+4,0],file[i+5,0],file[i+6,0],file[i+7,0]]
				for j in range(1,(column_lenght(elements,6))):
					if merger(temp_element) == grid.loc[j,'6_fold_index']:
						grid.loc[j,'6_fold'] += 1
		grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')]=(grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')])/Nframe
		plt.figure(x)
		sns.heatmap(grid.loc[:,('3_fold','4_fold','Tetrahedral','5_fold','6_fold')],
			annot=(grid.loc[:,('3_fold_index','4_fold_index','Tetrahedral_index','5_fold_index','6_fold_index')]),
			linewidths=.5,cmap="viridis",fmt='',square=False, yticklabels=False)
		plt.title("Coordination Heatmap with %s as center atom" % x)
	#--------------------------------------------------------------------------------------------------------
	plt.show()