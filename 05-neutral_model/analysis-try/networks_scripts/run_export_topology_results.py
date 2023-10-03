'''
	This script is intended to export all plots with topologycal measures
	and computations. 

'''


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def analyze_fname(a):
	'''
	Function that looks in file_name the values of
	temperature and step
	
	Returns:  
		- temperature (Celsius degrees)
		- step
	'''
	if a[0] == 'c':
		step = int(a.split('-')[1][1:])
		t=43
	elif a[0] == '3':
		step = int(a.split('-')[2])
		t=30
	return t, step

def redact_fname(temp, step):
	'''
        Function that gives the filename, having the temperature and step desired
    '''
	if step<10:
		step_string = '0' + str(step)
	elif step>=10:
		step_string = str(step)
	
	if temp == 30:
		x='30-1-'+step_string+'-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed_Nrem_rem'
		return x
	elif temp == 43:
		x='c43-p'+step_string+'-3-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed_Nrem_rem'
		return x
	else:
		print(f'No temperature {temp}')
		return None


def calculateParametersLSR(x, y):
	"""
	Function that calculates w* and b* for a least square regression and the standard error for b*
	Inputs:alpha = 0.5, color = "green", label = "$P(k)$")
plt.plot(k_, pkBin_fit, linestyle="dashed", color = "maroon", 
		label = "fit
		x: array containing the x-values in log-scale
		y: array containing the y-values in log-scale
	Outputs:
		b: float containing the intercept term of the regression line
		w: float containing the slope of the regression line
		SEw: float containing the standard error of b*
	"""

	meanTarget = np.mean(y)
	meanFeature = np.mean(x)

	centeredTarget = y - meanTarget
	centeredFeature = x - meanFeature

	w = (centeredFeature @ centeredTarget)/(centeredFeature @ centeredFeature)

	b = meanTarget - w *  meanFeature

	# Standard error
	yHat = b + w * x
	n = x.shape[0]

	SEw = np.sqrt((1/(n - 2)) * ((np.sum(np.power((y - yHat), 2)))/(np.sum(np.power((x - meanFeature), 2)))))



	return b, w, SEw

def histogramDegrees(degrees, minDegreeInterval = 0, maxDegreeInterval = 10, base = 2.0):
	"""
	Function that performs a logarithmic binning for the degrees in a network 
	Inputs:
		degrees: array of ints containing the degrees of a network 
		minDegreeInterval: int indicating the minimum degree to create the binning intervals
			Beware! this depends on the base because minDegree = base^minDegreeInterval
		maxDegreeInterval: int indicating the maximum degree to create the binning intervals
			Beware! this depends on the base because maxDegree = base^maxDegreeInterval 
		base: int/float indicating the base for the logarithmic binning
	Outputs:
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width
		binsUnique: array of ints (since intervals are discrete) containing the intervals of all bins 
		width: array of ints/floats containing the width of each interval 
	"""

	# Create array of bins
	bins = np.logspace(minDegreeInterval, maxDegreeInterval, base = base, dtype = int)

	# Get the unique bins
	binsUnique = np.unique(bins)
	# Get the widths of each interval
	widths = (binsUnique[1:] - binsUnique[:-1])
	# Get the histogram
	histogram = np.histogram(degrees, bins = binsUnique)

	# Normalize bins by width 
	histogram_temp = histogram[0]/widths

	# Normalize again by the number of nodes in the network, so we have a probability distribution
	histogram_norm = histogram_temp/degrees.shape[0]

	return histogram_norm, binsUnique, widths

def alphaNewmanClauset(degrees, kmin):
	"""
	Function that calculates the exponent coefficient (alpha) for the degree distribution of a network
	Inputs: 
		degrees: array of ints containing the degrees of a network
		kmin: int indicating the minimum degree for which the power law holds
	Outputs:
		alpha: float containing the exponent of the power law
			alpha = 1 + N * (sum_{i}^{N}log((k_i)/(k_min - 0.5)))^-1
		sigma: float containing the statistical error of the estimation of alpha 
			sigma = (alpha - 1)/(srqt(N))
	"""

	# Get the number of nodes with a degree greater than or equal to kmin
	#N = degrees[degrees >= kmin].shape[0]
	degrees = np.array(degrees)
	k_min = degrees[degrees >= kmin]
	N = k_min.shape[0]

	# Perform the sum
	summ = np.sum(np.log((k_min)/(kmin - 0.5)))

	# Get alpha
	alpha = 1 + N * (1/summ)
	# Get sigma
	sigma = (alpha - 1)/np.sqrt(N)

	return alpha, sigma

def histogramClustCoeffKnnAlternative(parameterDict, bins, widths):
	"""
	Function that performs a logarithmic binning over the nodes of a network by
	building a histogram and normalizing each bin by its width. Beware! This functions 
	only works for C(k) and Knn(k).
	Inputs:	
		parameterDict: dictionary whose keys are degrees and values are lists of C(k) or Knn(k) 
				k: [C_k_1, C_k_2, C_k_3]  -> k: [0, 0.33, 0.5]  ...
		bins: array of ints containing the intervals for all bins. The intervals are discrete.
		widts: array of ints: containing the width of each interval. 
	Outputs: 
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width.
	"""

	# Declare some variables
	#print('alternative histogram')
	histogram = []

	# Iterate over the bins
	for i in range(bins.shape[0] - 1):
		# The intervals
		minDegree = bins[i]
		# -1 since the right extreme of the interval is not taken into account and degrees are ints
		maxDegree = bins[i + 1] - 1

		#print("minDegree", minDegree, "maxDegree", maxDegree)

		# Count all nodes that fall in the interval
		if minDegree == maxDegree:
			# Sum the parameter of all the point in the current interval
			suma = 0
			
			#print("\tsumaa", suma)
			for paramPoint in parameterDict[minDegree]:
				suma += paramPoint
			# Count the number of nodes at the current interval
			n = len(parameterDict[minDegree])
			histogram.append(suma/n)

		else:
			count = 0
			suma = 0
			# In case the interval comprises more than one degree
			#print("interval comprises more than one degree\n")
			for j in range(minDegree, maxDegree + 1):
				#  Check if the degree exists
				if j in parameterDict.keys():
					#print("\tCheck the degree in interval", i)
					# Count the number of nodes at the current interval
					count += len(parameterDict[j])
					#print("counttt", count, "degree", i)
					# Sum the parameter of all the point in the current interval
					#suma += parameterDict[j]
					for paramPoint in parameterDict[j]:
						suma += paramPoint


			# In case no points fall in the interval
			#print("\tsumaa", suma)
			if suma == 0:
				histogram.append(0)
			else:
				histogram.append(suma/count)

	# Normalize histogram
	histogram_norm = np.array(histogram)#/widths

	return histogram_norm

def histogramClustCoeffKnn(parameterDict, bins, widths):
	"""
	Function that performs a logarithmic binning over the nodes of a network by
	building a histogram and normalizing each bin by its width. Beware! This functions 
	only works for C(k) and Knn(k).
	Inputs:	
		parameterDict: dictionary whose keys are degrees and values are lists that contains 
		nodes with a parameter associated to a degree.
		bins: array of ints containing the intervals for all bins. The intervals are discrete.
		widts: array of ints: containing the width of each interval. 
	Outputs: 
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width.
	"""

	# Declare some variables
	histogram = []

	print('\ncreating histogram!')

	# Iterate over the bins
	for i in range(bins.shape[0] - 1):
		# The intervals
		minDegree = bins[i]
		# -1 since the right extreme of the interval is not taken into account and degrees are ints
		maxDegree = bins[i + 1] - 1

		print("minDegree", minDegree, "maxDegree", maxDegree)

		# Count all nodes that fall in the interval
		if minDegree == maxDegree:
			# Sum the parameter of all the point in the current interval
			suma = 0 
			for paramPoint in parameterDict[minDegree]:
				suma += paramPoint
			# Count the number of nodes at the current interval
			n = len(parameterDict[minDegree])
			histogram.append(suma/n)

		else:
			count = 0
			suma = 0
			# In case the interval comprises more than one degree
			print("interval comprises more than one degree\n")
			for j in range(minDegree, maxDegree + 1):
				#  Check if the degree exists
				if j in parameterDict.keys():
					print("\tCheck the degree in interval", i)
					# Count the number of nodes at the current interval
					count += len(parameterDict[j])
					#print("counttt", count, "degree", i)
					# Sum the parameter of all the point in the current interval
					for paramPoint in parameterDict[j]:
						suma += paramPoint

			print("\tcount outside", count)
			print("\tsum outside", suma)
			# In case no points fall in the interval
			if count == 0:
				print("\tno points fall in the interval!!!")
				histogram.append(0)
			else:
				histogram.append(suma/count)

		print("iii", i)

	# Normalize histogram
	histogram_norm = np.array(histogram)/widths

	print('shapes', len(histogram), len(widths))

	return histogram_norm

#top_res_path = './Topology_results/'
top_res_path = './Topology_results/Progressive_GN/'
temps = [30, 43]
cmap = {30:'blue', 43:'red'}
symb_dict = {2:"o" , 4: "^", 10: "s", 30: "p", 60: "H"}
greek_letterz=[chr(code) for code in range(945,970)]


#Results of
#GN : number of nodes, number of edges, and number of nodes &edges in Giant component
#Repeat the same for each step and temp
#Save a file like this:
#GN, -, -, fraction,-, -, fraction
#t_step, -, -, fraction,-, -, fraction


#Activate this to export the evolution of mean degree and std.
if 1==1:
	print('#############Starting plot of mean degrees#####################')
	############ Degree #############
	#Read the graph
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4, 10, 30, 60]

	#mean k vs step for both temps in a plot
	#std-dev of k vs step for both temps in a plot
	results = {t:{} for t in temps}
	
	for thisFileName in allFileNames:
		try:
			t, step = analyze_fname(thisFileName.replace('edges_', ''))
		except:
			t, step = analyze_fname(thisFileName.replace(' edges_', ''))

		df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
		results[t][step] = {'mean': df['degree'].mean(), 'std': df['degree'].std()}
		
	fig, axs = plt.subplots(1,2, figsize=(10, 5))
	for t in temps:
		L = list(results[t].keys())
		L.sort()

		x = np.array(L)
		y = []
		y_err = []
		for step in L:
			y.append(results[t][step]['mean'])
			y_err.append(results[t][step]['std'])
		y = np.array(y)
		y_err = np.array(y_err)

		#b,w,E = calculateParametersLSR(np.log2(x), np.log2(y))
		#y_aster = 2**b * x**w 


		axs[0].plot(x,y, 'o-', c=cmap[t],  label=f'{t}ºC')
		#axs[0].plot(x,y_aster, '--', c='black', label=f'fiting {round(w,3)}({int(E*1000)})')
		axs[0].set_xlabel('step')
		axs[0].set_ylabel('<k>')
		axs[0].legend(loc='best')
		axs[1].plot(x,y_err, 'o-', c=cmap[t],  label=f'{t}ºC')
		axs[1].set_xlabel('k')
		axs[1].set_ylabel(f"{greek_letterz[18]}(k)")
		axs[1].legend(loc='best')
		#axs[0].set_yscale('log')
		#axs[0].set_xscale('log')
		#axs[1].set_yscale('log')
		#axs[1].set_xscale('log')
	#Add GN as gridline
		
	fig.tight_layout()
	plt.savefig(f'{top_res_path}/Figures/k.svg', dpi=300, format='svg')

 

	fig, axs = plt.subplots(1,2, figsize=(10,5))

	#plot of p(k) vs k for both temps and steps=[6, 10,30,60] & GN & ancestor (step=2)
	for thisFileName in allFileNames:
		p_k = {} 
		try:
			t, step = analyze_fname(thisFileName.replace('edges_', ''))
		except:
			t, step = analyze_fname(thisFileName.replace(' edges_', ''))
		
		if t==30:
			a = 0
		else:
			a=1
		if step in desired_steps:
			df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			k_list = df['degree'].tolist()
			for k in k_list:
				try:
					p_k[k] += 1
				except:
					p_k[k] = 1
			normalization =np.sum(list(p_k.values()))

			x = np.array(list(p_k.keys()))

			y = np.array([p_k[k]/normalization for k in x])
			axs[a].plot(x,y, marker=symb_dict[step] , linestyle='none', alpha=0.7, label=f"{step}")
			axs[a].set_xscale('log')
			axs[a].set_yscale('log')
			
			

			axs[a].set_xlabel('k')
			axs[a].set_ylabel('p(k)')
			axs[a].legend(loc='best')
			axs[a].set_title(f'{t}ºC')
	#Add ancestor and GN
	for t in temps:
		p_k = {} 
		p_k_GN = {} 
		step=2
		fname = redact_fname(temp=t, step=2) + '.csv'
		#df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
		
		if t==30:
			a = 0
		else:
			a = 1
		if 1==1:
			print('Ancestor', fname)
			df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
			df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
			k_list = df['degree'].tolist()
			k_list_GN = df_GN['degree'].tolist()
			for k in k_list:
				try:
					p_k[k] +=1    
				except:
					p_k[k] = 1    

			normalization =np.sum(list(p_k.values()))
			
			for k in k_list_GN:
				try:
					p_k_GN[k] +=1
				except:
					p_k_GN[k] = 1
			normalization_GN =np.sum(list(p_k_GN.values()))

			#for k in k_list:
			#    p_k[k] = p_k[k]/normalization
			#print(p_k)
			x = np.array(list(p_k.keys()))
			x_GN = np.array(list(p_k_GN.keys()))

			y = np.array([p_k[k]/normalization for k in x])
			#print(y.sum())
			y_GN = np.array([p_k_GN[k]/normalization_GN for k in x_GN])
			#print(y_GN.sum())
			axs[a].plot(x,y, marker=symb_dict[step], fillstyle='none', linestyle='none', c='black', alpha=0.8,label=f"Ancestor")
			axs[a].plot(x_GN,y_GN, 'o', c='gray', linestyle='none', alpha=0.3,label=f"Grand Network")
			axs[a].set_xscale('log')
			axs[a].set_yscale('log')
			
			

			axs[a].set_xlabel('k')
			axs[a].set_ylabel('p(k)')
			axs[a].legend(loc='best')
		

	fig.tight_layout()
	plt.savefig(f'{top_res_path}Figures/pk.svg', dpi=300, format='svg')   


#Activate this to export the cummulative p(k)
if 0==1:
	print('#############Starting plot of cummulative P(k) #####################')
	for t_desired in temps:
		allFileNames = os.listdir(f'{top_res_path}Neighbors/')
		allFileNames.sort()
		desired_steps = [4,10,30,60]

		#Initiate plot
		plt.figure()

		point_sizes = 30
		for thisFileName in allFileNames:
			p_k = {}
			t, step = analyze_fname(thisFileName.replace('edges_', ''))
			
			if step in desired_steps and t==t_desired:
				df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
				degrees = np.array(df['degree'].tolist())
				values, counts = np.unique(degrees, return_counts=True)
				counts = counts/counts.sum()
				counts_cum = counts.copy()
				for i in range(1,counts_cum.shape[0]):
					counts_cum[i] = counts_cum[i-1] + counts[i]
				counts_cum = 1.0-1.0*counts_cum
				plt.scatter(values, counts_cum, alpha = 0.5, marker=symb_dict[step],label=f'{step}')
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('k')
		plt.ylabel('1-Cummulative p(k)')
		plt.legend(loc='best')
		plt.ylim(10**-7,1.5)
		plt.savefig(f'{top_res_path}Figures/pk_cummulative_{t_desired}.svg', dpi=300, format='svg')



#Activate this to compute selected steps, ancestor and GN, and fit these last two
if 0==1:
	print('#############Starting plot of selected steps#####################')
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]
	fit_min = 3
	fit_max = 25

	#Initiate plot
	plt.figure()
	point_sizes = 30
	for thisFileName in allFileNames:
		#p_k = {}
		t, step = analyze_fname(thisFileName.replace('edges_', ''))
        

		df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
		degrees = np.array(df['degree'].tolist())

		if step in desired_steps and t==43:               #Selected steps by Henry
			print("Processing "+thisFileName);
			# Logarithmic binning
			histDegrees_norm, bins, widths = histogramDegrees(degrees)
			k_ = bins[:-1]
			pkBin_ = histDegrees_norm
			# Least squares regression
			k_scaled = np.log10(k_[pkBin_ != 0])
			pkBin_scaled = np.log10(pkBin_[pkBin_ != 0])
			
			plt.xscale('log', nonpositive='clip')
			plt.yscale('log', nonpositive='clip')
			plt.scatter(k_, pkBin_, alpha = 0.5, marker=symb_dict[step], label = f"{step}")
			
			

			
	
	
	
	#Add ancestor
	step=2
	fname = redact_fname(temp=t, step=2) + '.csv'
	df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	degrees = np.array(df['degree'].tolist())
	degrees = degrees[degrees >=2]
	# Logarithmic binning
	histDegrees_norm, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]
	pkBin_ = histDegrees_norm
	# Least squares regression
	k_scaled = np.log10(k_[pkBin_ != 0])
	pkBin_scaled = np.log10(pkBin_[pkBin_ != 0])
	b, w, SEw = calculateParametersLSR(k_scaled[1:], pkBin_scaled[1:])
	pkBin_fit = 10 ** (b) * k_** w	
	plt.plot(k_, pkBin_fit, linestyle="--", color = "grey", label = "Ancestor alpha: " + str(round(w, 2)) + "$\\pm: $" + str(round(SEw, 3)))
	plt.xscale('log', nonpositive='clip')
	plt.yscale('log', nonpositive='clip')
	plt.scatter(k_, pkBin_, alpha = 0.5, marker='o', facecolors='none', edgecolors='black', label = f"Ancestor")
	

	

	#Add GrandNetwork
	p_k_GN = {}
	df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	degrees = np.array(df_GN['degree'].tolist())
	degrees = degrees[degrees >=3]
	# Logarithmic binning
	histDegrees_norm, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]
	pkBin_ = histDegrees_norm
	# Least squares regression
	k_scaled = np.log10(k_[pkBin_ != 0])
	pkBin_scaled = np.log10(pkBin_[pkBin_ != 0])
	b, w, SEw = calculateParametersLSR(k_scaled[1:], pkBin_scaled[1:])
	pkBin_fit = 10 ** (b) * k_** w	
	plt.plot(k_, pkBin_fit, linestyle="-", color = "gray", label = "Grand Network alpha: " + str(round(w, 2)) + "$\\pm: $" + str(round(SEw, 3)))
	plt.xscale('log', nonpositive='clip')
	plt.yscale('log', nonpositive='clip')
	plt.scatter(k_, pkBin_, alpha = 0.5, marker='o', c='gray', edgecolor='black', label = f"Grand Network")

	plt.legend()
	plt.xlabel('k')
	plt.ylabel('p(k)')

	plt.tight_layout()
	plt.savefig(f'{top_res_path}Figures/pk_binned.svg', dpi=300, format='svg')
	



#Activate this to compute alphas in all steps
if 0==1:
	print('#############Starting plot of alphas#####################')
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	k_min = 2
	k_min_GN = 3
	results_pareto = {t:{} for t in temps}

	plt.figure()
	avoid_steps = [1,2,4,10,60]
	#avoid_steps = []
	alphas = {t: {} for t in temps}
	alphas_err = {t: {} for t in temps}
	for thisFileName in allFileNames:
		
		t, step = analyze_fname(thisFileName.replace('edges_', ''))
        

		df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
		results_pareto[t][step] = {'mean': df[df['degree']>=k_min]['degree'].mean(), 'std': df[df['degree']>=2]['degree'].std()}
		print("Processing "+thisFileName); 



		if step:               #Selected steps by Henry
			k_list = df['degree'].tolist()
			alpha, sigma = alphaNewmanClauset(k_list, k_min)
			alphas[t][step]=alpha
			alphas_err[t][step]=sigma
	
	
	
	#Add GrandNetwork
	df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	k_list_GN = df_GN['degree'].tolist()
	k_mean_GN = np.array(k_list_GN)[np.array(k_list_GN)>=k_min_GN].mean()
	alpha_GN, sigma_GN = alphaNewmanClauset(k_list_GN, k_min_GN)

	for t in temps:
		x = list(alphas[t].keys())
		x.sort()
		y = np.array([-1.0*alphas[t][s] for s in x])
		y_err = np.array([alphas_err[t][s] for s in x])
		plt.plot(x,y, '-o' ,c=cmap[t], label=f't={t}')
		plt.fill_between(x, y-y_err, y+y_err, alpha=0.5)
	plt.plot(x, [-1.0*alpha_GN for element in x], '--', c='black', label='Grand Network')
	
	#Add ancestor
	fname = redact_fname(temp=43, step=2) + '.csv'
	df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	
	k_list = df['degree'].tolist()
	alpha, sigma = alphaNewmanClauset(k_list, k_min)
	plt.plot(x, [-1.0*alpha for element in x], 'o-', c='black', label='Ancestor')


	plt.legend(loc='best')
	plt.ylabel(f'{chr(945)}')
	plt.xlabel('step')
	plt.savefig(f'{top_res_path}Figures/alpha.svg', dpi=300, format='svg')

	plt.figure()

	for t in temps:
		L = list(results[t].keys())
		L.sort()
		y = [alphas[t][step] for step in L]
		y_err = [alphas_err[t][step] for step in L]
		x = [results_pareto[t][step]['mean'] for step in L]
		x_err = [results_pareto[t][step]['std'] for step in L]
		#b,w,er = calculateParametersLSR(x, y)
		plt.errorbar(x,y, yerr=y_err, c=cmap[t], linestyle='', marker='o', label=f'{t}')
		#x_aster = np.linspace(min(x), max(x), 100)
		#y_aster = w*x_aster + b
		#plt.plot(x_aster,y_aster, '--', c=cmap[t], label=f'Fitting slope={round(w,2)} intercept={round(b,2)}')
	#y_aster_aster = (x_aster/2)/((x_aster/2)-1)+1 - 0.4
	#plt.plot(x_aster,-1.0*y_aster_aster, '--', c='green', label=f'Sample')
	plt.legend()
	plt.xlabel('<k>=2> ')
	plt.ylabel(f'{chr(945)}')
	plt.savefig(f'{top_res_path}Figures/alpha_k_dependence.svg', dpi=300, format='svg')


############# Assortativity #############
if 0==1:
	#plot of Knn vs. k for both temps and steps=[6,10,30,60] & GN & ancestor (step=2). Add fitting in range k=[7, 200] for GN & ancestor
	#labels K_nn(k) vs. k in log-log 

	print('#############Starting plot of Knn vs. k#####################')
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]

	min_bin_fit = 7
	max_bin_fit = 30

	#Initiate plot
	plt.figure()
	
	point_sizes = 30
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step in desired_steps and t==43:               #Selected steps by Henry
			print("Processing "+thisFileName);

			df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			degrees = np.array(df['degree'].tolist())
			Knn = np.array(df['average_neighbor_degree'].tolist())
			p_dict = {k: [] for k in set(list(degrees))}

			for k,nn in zip(degrees, Knn):
				p_dict[k].append(nn)
			
			#Obtain bins and widths
			pkBin_, bins, widths = histogramDegrees(degrees)
			k_ = bins[:-1]

			#Obtain the re-binning
			KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

			

			
			plt.xscale('log', nonpositive='clip')
			plt.yscale('log', nonpositive='clip')
			plt.scatter(k_, KNN_bin, alpha = 0.5, marker=symb_dict[step], label = f"{step}")
			#print(f'Fitting: slope={pKNNw}+-{pKNNESw}')
	

	#Add ancestor
	step=2
	fname = redact_fname(temp=t, step=2) + '.csv'
	df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	degrees = np.array(df['degree'].tolist())
	Knn = np.array(df['average_neighbor_degree'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw, pKNNESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

	# Get regression lines and add to list
	KNN_bin_fit = 10 ** (pKNNb) * k_** pKNNw

	plt.plot(k_, KNN_bin_fit, linestyle="--", color = "grey", label = "Ancestor alpha: " + str(round(pKNNw, 2)) + "$\\pm: $" + str(round(pKNNESw, 3)))
	plt.scatter(k_, KNN_bin, alpha = 0.5,marker='o', facecolors='none', edgecolors='black', label = f"Ancestor")
	

	#Add GrandNetwork
	p_k_GN = {}
	df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	degrees = np.array(df_GN['degree'].tolist())
	Knn = np.array(df_GN['average_neighbor_degree'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw, pKNNESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

	# Get regression lines and add to list
	KNN_bin_fit = 10 ** (pKNNb) * k_** pKNNw

	plt.plot(k_, KNN_bin_fit, linestyle="-", color = "gray", label = "Grand Network alpha: " + str(round(pKNNw, 2)) + "$\\pm: $" + str(round(pKNNESw, 3)))

	plt.scatter(k_,KNN_bin, alpha = 0.5, marker='o', c='gray', edgecolor='black', label = f"Grand Network")

	plt.legend()
	plt.xlabel('k')
	plt.ylabel('Knn(k)')

	plt.tight_layout()
	
	plt.savefig(f'{top_res_path}Figures/Knn_k.svg', dpi=300, format='svg')


#Activate this to compute the exponent of power law region in Knn vs. k plots
if 0==1:
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]

	min_bin_fit = 7
	max_bin_fit = 30

	slopes = {t: {} for t in temps}
	slopes_err = {t: {} for t in temps}

	#Initiate plot
	plt.figure()
	
	point_sizes = 30
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step and t:           
			print("Processing "+thisFileName);

			df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			degrees = np.array(df['degree'].tolist())
			Knn = np.array(df['average_neighbor_degree'].tolist())
			p_dict = {k: [] for k in set(list(degrees))}

			for k,nn in zip(degrees, Knn):
				p_dict[k].append(nn)
			
			#Obtain bins and widths
			pkBin_, bins, widths = histogramDegrees(degrees)
			k_ = bins[:-1]

			#Obtain the re-binning
			KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

			#Fitting
			k_scaled = np.log10(k_[KNN_bin != 0])
			KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
			pKNNb, pKNNw, pKNNESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

			slopes[t][step] = pKNNw
			slopes_err[t][step] = pKNNESw	
	
	
	#Add ancestor
	step=2
	fname = redact_fname(temp=t, step=2) + '.csv'
	df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	degrees = np.array(df['degree'].tolist())
	Knn = np.array(df['average_neighbor_degree'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}
	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw_anc, pKNNESw_anc = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

	#Add GrandNetwork
	p_k_GN = {}
	df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	degrees = np.array(df_GN['degree'].tolist())
	Knn = np.array(df_GN['average_neighbor_degree'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw_GN, pKNNESw_GN = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])


	for t in temps:
		L = list(slopes[t].keys())
		L.sort()
		y = [slopes[t][step] for step in L]
		y_err = [slopes_err[t][step] for step in L]
		x = np.array(L)
		plt.errorbar(x,y, yerr=y_err, c=cmap[t], linestyle='', marker='o', label=f'{t}')
		
	plt.plot(x, [pKNNw_anc for i in x], '--', label='Ancestor')
	plt.plot(x, [pKNNw_GN for i in x], '-', label='Grand Network')

	plt.legend()
	plt.xlabel('step ')
	plt.ylabel(f'{chr(945)}')
	plt.savefig(f'{top_res_path}Figures/Knn_k_slopes.svg', dpi=300, format='svg')


############# Clusstering #############
if 0==1:
	#plot of C(k) vs. k for both temps and steps=[6,10,30,60] & GN & ancestor (step=2). Add fitting in range k=[7, 200] for GN & ancestor
	#labels C(k) vs. k in log-log 

	print('#############Starting plot of C(k) vs. k#####################')
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]
	#desired_steps = [12, 14, 16, 18, 20]
	temp = 43

	min_bin_fit = 4
	max_bin_fit = 21

	#Initiate plot
	plt.figure()
	
	point_sizes = 30
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step in desired_steps and t==temp:               #Selected steps by Henry
			print("Processing "+thisFileName);

			#df_1 = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			df = pd.read_csv(f'{top_res_path}Clustering Coefficients/{thisFileName}')
			#df = pd.merge(df_1, df_2, on='node')
			
			degrees = np.array(df['degree'].tolist())
			
			Ck = np.array(df['LCC'].tolist())
			p_dict = {k: [] for k in set(list(degrees))}

			for k,nn in zip(degrees, Ck):
				p_dict[k].append(nn)
			
			#Obtain bins and widths
			pkBin_, bins, widths = histogramDegrees(degrees)
			k_ = bins[:-1]

			#Obtain the re-binning
			Ck_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)
					
			plt.xscale('log', nonpositive='clip')
			plt.yscale('log', nonpositive='clip')
			try: plt.scatter(k_, Ck_bin, alpha = 0.5, marker=symb_dict[step], label = f"{step}") #
			except: plt.scatter(k_, Ck_bin, alpha = 0.5, label = f"{step}")
			#print(f'Fitting: slope={pKNNw}+-{pKNNESw}')
	

	#Add ancestor
	step=2
	fname = redact_fname(temp=temp, step=2) + '.csv'
	#df_1 = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	df = pd.read_csv(f'{top_res_path}Clustering Coefficients/edges_{fname}')
	#df = pd.merge(df_1, df_2, on='node')
	degrees = np.array(df['degree'].tolist())
	Ck = np.array(df['LCC'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Ck):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#print('esto ', [f'{i} {k}' for i,k in enumerate(k_)])

	#Obtain the re-binning
	Ck_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[Ck_bin != 0])
	Ck_bin_scaled = np.log10(Ck_bin[Ck_bin != 0])
	pCkb, pCkw, pCkESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], Ck_bin_scaled[min_bin_fit:max_bin_fit])

	# Get regression lines and add to list
	Ck_bin_fit = 10 ** (pCkb) * k_** pCkw

	plt.plot(k_, Ck_bin_fit, linestyle="--", color = "grey", label = "Ancestor alpha: " + str(round(pCkw, 2)) + "$\\pm: $" + str(round(pCkESw, 3)))
	plt.scatter(k_, Ck_bin, alpha = 0.5,marker='o', facecolors='none', edgecolors='black', label = f"Ancestor")
	

	#Add GrandNetwork
	p_k_GN = {}
	#df_GN_1 = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	df = pd.read_csv(f'{top_res_path}clustering_GrandNetwork.csv')
	#df = pd.merge(df_GN_1, df_GN_2, on='node')
	degrees = np.array(df['degree'].tolist())
	Ck = np.array(df['LCC'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Ck):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	Ck_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[Ck_bin != 0])
	Ck_bin_scaled = np.log10(Ck_bin[Ck_bin != 0])
	pCkb, pCkw, pCkESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], Ck_bin_scaled[min_bin_fit:max_bin_fit])

	# Get regression lines and add to list
	Ck_bin_fit = 10 ** (pCkb) * k_** pCkw

	plt.plot(k_, Ck_bin_fit, linestyle="-", color = "gray", label = "Grand Network alpha: " + str(round(pCkw, 2)) + "$\\pm: $" + str(round(pCkESw, 3)))

	plt.scatter(k_,Ck_bin, alpha = 0.5, marker='o', c='gray', edgecolor='black', label = f"Grand Network")

	plt.legend()
	plt.xlabel('k')
	plt.ylabel('C(k)')

	plt.tight_layout()
	
	plt.savefig(f'{top_res_path}Figures/Ck_k_{temp}.svg', dpi=300, format='svg')


#Activate this to compute the exponent of power law region in C(k) vs. k plots
if 0==1:
	allFileNames = os.listdir(f'{top_res_path}Neighbors/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]

	min_bin_fit = 4
	max_bin_fit = 21

	slopes = {t: {} for t in temps}
	slopes_err = {t: {} for t in temps}

	#Initiate plot
	plt.figure()
	
	point_sizes = 30
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step and t:           
			print("Processing "+thisFileName);

			#df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			
			df = pd.read_csv(f'{top_res_path}Clustering Coefficients/{thisFileName}')
			degrees = np.array(df['degree'].tolist())
			Knn = np.array(df['LCC'].tolist())
			p_dict = {k: [] for k in set(list(degrees))}

			for k,nn in zip(degrees, Knn):
				p_dict[k].append(nn)
			
			#Obtain bins and widths
			pkBin_, bins, widths = histogramDegrees(degrees)
			k_ = bins[:-1]

			#Obtain the re-binning
			KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

			
			k_scaled = np.log10(k_[KNN_bin != 0])
			KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
			pKNNb, pKNNw, pKNNESw = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

			slopes[t][step] = pKNNw
			slopes_err[t][step] = pKNNESw	
	
	
	#Add ancestor
	step=2
	fname = redact_fname(temp=t, step=2) + '.csv'
	#df = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	
	df = pd.read_csv(f'{top_res_path}Clustering Coefficients/edges_{fname}')
	degrees = np.array(df['degree'].tolist())
	Knn = np.array(df['LCC'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}
	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw_anc, pKNNESw_anc = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])

	#Add GrandNetwork
	p_k_GN = {}
	#df_GN = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	
	df_GN = pd.read_csv(f'{top_res_path}clustering_GrandNetwork.csv')
	degrees = np.array(df_GN['degree'].tolist())
	Knn = np.array(df_GN['LCC'].tolist())

	p_dict = {k: [] for k in set(list(degrees))}

	for k,nn in zip(degrees, Knn):
		p_dict[k].append(nn)
			
	#Obtain bins and widths
	pkBin_, bins, widths = histogramDegrees(degrees)
	k_ = bins[:-1]

	#Obtain the re-binning
	KNN_bin = histogramClustCoeffKnnAlternative(p_dict, bins, widths)

	#Fitting
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	pKNNb, pKNNw_GN, pKNNESw_GN = calculateParametersLSR(k_scaled[min_bin_fit:max_bin_fit], KNN_bin_scaled[min_bin_fit:max_bin_fit])


	for t in temps:
		L = list(slopes[t].keys())
		L.sort()
		y = [slopes[t][step] for step in L]
		y_err = [slopes_err[t][step] for step in L]
		x = np.array(L)
		plt.errorbar(x,y, yerr=y_err, c=cmap[t], linestyle='', marker='o', label=f'{t}')
		
	plt.plot(x, [pKNNw_anc for i in x], '--', label='Ancestor')
	plt.plot(x, [pKNNw_GN for i in x], '-', label='Grand Network')

	plt.legend()
	plt.xlabel('step ')
	plt.ylabel(f'{chr(945)}')
	plt.savefig(f'{top_res_path}Figures/Ck_k_slopes.svg', dpi=300, format='svg')


######## Centrality ##############

############ Betweenness centrality ##############

if 0==1:
	#plot of C(k) vs. k for both temps and steps=[6,10,30,60] & GN & ancestor (step=2). Add fitting in range k=[7, 200] for GN & ancestor
	#labels C(k) vs. k in log-log 

	print('#############Starting plot of B(i) vs. k#####################')
	allFileNames = os.listdir(f'{top_res_path}Betweenness_centrality/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]
	#desired_steps = [12, 14, 16, 18, 20]
	temp = 43

	min_bin_fit = 4
	max_bin_fit = 21

	#Initiate plot
	plt.figure()
	
	point_sizes = 20
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step in desired_steps and t==temp:               #Selected steps by Henry
			print("Processing "+thisFileName);

			df_1 = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			df_2 = pd.read_csv(f'{top_res_path}Betweenness_centrality/{thisFileName}')
			df = pd.merge(df_1, df_2, on='node')
			
			try: plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5, marker=symb_dict[step], label = f"{step}") #
			except: plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5, label = f"{step}")

	
	
	#Add ancestor
	step=2
	fname = redact_fname(temp=temp, step=2) + '.csv'
	df_1 = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	df_2 = pd.read_csv(f'{top_res_path}Betweenness_centrality/edges_{fname}')
	df = pd.merge(df_1, df_2, on='node')
	plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5,marker='o', facecolors='none', edgecolors='black', label = f"Ancestor")
	plt.legend()
	plt.xlabel('k')
	plt.ylabel('B(i)')
	plt.yscale('log')
	plt.xscale('log')

	plt.tight_layout()


	#Add GrandNetwork    This is commented because we don't have the measures
	#p_k_GN = {}
	#df_1 = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	#df_2 = pd.read_csv(f'{top_res_path}betweenness_GrandNetwork.csv')
	#df = pd.merge(df_1, df_2, on='node')
	#plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5, marker='o', c='gray', edgecolor='black', label = f"Grand Network")
	
	plt.savefig(f'{top_res_path}Figures/Bi_k_{temp}.png', dpi=300, format='png')

########### Eigenvector centrality ##############
if 0==1:
	#plot of C(k) vs. k for both temps and steps=[6,10,30,60] & GN & ancestor (step=2). Add fitting in range k=[7, 200] for GN & ancestor
	#labels C(k) vs. k in log-log 

	print('#############Starting plot of E	(i) vs. k#####################')
	allFileNames = os.listdir(f'{top_res_path}Eigenvector_centrality/')
	allFileNames.sort()
	desired_steps = [4,10,30,60]
	#desired_steps = [12, 14, 16, 18, 20]
	temp = 43

	min_bin_fit = 4
	max_bin_fit = 21

	
	
	point_sizes = 20
	for thisFileName in allFileNames:
		t, step = analyze_fname(thisFileName.replace('edges_', ''))

		if step in desired_steps and t==temp:               #Selected steps by Henry
			print("Processing "+thisFileName);
			#Initiate plot
			plt.figure()

			df_1 = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
			df_2 = pd.read_csv(f'{top_res_path}Eigenvector_centrality/{thisFileName}')
			df = pd.merge(df_1, df_2, on='node')
			
			try: plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5, marker=symb_dict[step]) #
			except: plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5)
			#plt.legend()
			plt.xlabel('k')
			plt.ylabel('v1(i)')
			plt.yscale('log')
			plt.xscale('log')
			plt.ylim(10**-10,1)

			plt.tight_layout()
			plt.savefig(f'{top_res_path}Figures/ei_k_{temp}_{step}.png', dpi=300, format='png')
			#print(f'saved {top_res_path}Figures/ei_k_{temp}_{step}.png')

	
	
	#Add ancestor
	plt.figure()
	step=2
	fname = redact_fname(temp=temp, step=2) + '.csv'
	df_1 = pd.read_csv(f'{top_res_path}Neighbors/edges_{fname}')
	df_2 = pd.read_csv(f'{top_res_path}Eigenvector_centrality/edges_{fname}')
	df = pd.merge(df_1, df_2, on='node')
	plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5,marker='o', facecolors='none', edgecolors='black', label = f"Ancestor")
	#plt.legend()
	plt.xlabel('k')
	plt.ylabel('v1(i)')
	plt.yscale('log')
	plt.xscale('log')
	plt.ylim(10**-10,1)
	plt.tight_layout()
	plt.savefig(f'{top_res_path}Figures/ei_k_{temp}_anc.png', dpi=300, format='png')


	#Add GrandNetwork
	p_k_GN = {}
	#df = pd.read_csv(f'{top_res_path}neighbors_GrandNetwork.csv')
	df = pd.read_csv(f'{top_res_path}eigenvector_GrandNetwork.csv')
	#df = pd.merge(df_1, df_2, on='node')
	plt.scatter(df['degree'], df['eigenvector_centrality'], alpha = 0.5, marker='o', c='gray', edgecolor='black', label = f"Grand Network")
	#plt.legend()
	plt.xlabel('k')
	plt.ylabel('v1(i)')
	plt.yscale('log')
	plt.xscale('log')
	plt.ylim(10**-10,1)
	plt.savefig(f'{top_res_path}Figures/ei_k_{temp}_GN.png', dpi=300, format='png')





#For each temp: plot a 3rowsx2columns image containing: step 2 (ancestor), step 6,10, 30, 60 , & GN
#labels v_1(i) vs. k in log-log 

#plot abundances vs. eigenvector centrality. Also add an y=x line
#labels = abundance(i) vs. v_1(i) log-log
#Do the previous like a 3rowsx2columns image containing: step 2 (ancestor), step 6,10, 30, 60 , & GN

#Plot using Cyoscape of
##Size of the nodes proportional to eigenvector centrality
#image 1: circular graph of t=30 and step=6
#image 2: circular graph of t=30 and step=60
#image 3: circular graph of t=43 and step=6
#image 4: circular graph of t=43 and step=60


#pending analysis
# frequency of mutations vs. Position (loci)
# GN loci rank (by frequency of mutations). lin-lin, linx-logy , log-log
# GN loci within codon, frequency of mutation



#number of phenotypes encoded by a genotype  PHI(gamma) understand better
#