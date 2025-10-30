#!/usr/bin/python

'''
This program is designed to construct a plot comparing the SFS for multiple 
potential numbers of linkage blocks per chromosome, allowing for us to verify
that varying the number of linkage blocks doesn't have a significant impact
on the trends we see with increasing Ud. 

Author: Micaila Marcelle
'''

# Imports all necessary libraries/packages
import math
import tskit
import io
import numpy as np
import matplotlib.pyplot as plt
import msprime
import pandas as pd
from scipy.signal import savgol_filter

print("Successfully imported necessary modules")

cur_Ud = 10
list_of_Ls = [500, 1000, 1500, 2000]

siteFrequencySpectrum_500 = []
siteFrequencySpectrum_1000 = []
siteFrequencySpectrum_1500 = []
siteFrequencySpectrum_2000 = []

# Used for finding an "averaged out" Ne, since this will allow for a single baseline
# to compare the distortions to
tree_sequence_list_overall = []

for link_num in list_of_Ls:
      
	# Reads in the necessary tskit tables for all generations of interest, adding
	# the results to the corresponding list for each type of table
	node_list = []
	edge_list = []
	site_list = []
	mut_list = []

	# First set reads in the tskit tables for generation 500000
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/nodetable_gen500000.txt") as file:
		node_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/edgetable_gen500000.txt") as file:
		edge_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/sitetable_gen500000.txt") as file:
		site_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/mutationtable_gen500000.txt") as file:
		mut_list.append(file.read())
		
	# Second reads in the tskit tables for generation 1000000
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/nodetable_gen1000000.txt") as file:
		node_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/edgetable_gen1000000.txt") as file:
		edge_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/sitetable_gen1000000.txt") as file:
		site_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/mutationtable_gen1000000.txt") as file:
		mut_list.append(file.read())
		
	# Third reads in tskit tables for generation 1500000
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/nodetable_gen1500000.txt") as file:
		node_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/edgetable_gen1500000.txt") as file:
		edge_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/sitetable_gen1500000.txt") as file:
		site_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/mutationtable_gen1500000.txt") as file:
		mut_list.append(file.read())
		
	# Finally, we read in tskit tables for generation 2000000
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/nodetable.txt") as file:
		node_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/edgetable.txt") as file:
		edge_list.append(file.read())  
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/sitetable.txt") as file:
		site_list.append(file.read())
	with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/mutationtable.txt") as file:
		mut_list.append(file.read())
	'''
	else:
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/nodetable_gen2000000.txt") as file:
			node_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/edgetable_gen2000000.txt") as file:
			edge_list.append(file.read())  
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/sitetable_gen2000000.txt") as file:
			site_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_0.0095_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_" + str(link_num) + "_seed_24_Sd_0.009487/mutationtable_gen2000000.txt") as file:
			mut_list.append(file.read())
	'''
		
	print("Successfully opened table files")
    
	# Uses these tables in order to construct the corresponding tree sequences
	# Note that the process of constructing tree sequences is repeated for each
	# checkpoint, leading to a number of tree sequences equal to the number of
	# checkpoints
	tree_sequence_list = []
	for i in range(len(node_list)):
		tree_sequence_list.append(tskit.load_text(io.StringIO(node_list[i]), 
								io.StringIO(edge_list[i]), 
								io.StringIO(site_list[i]), 
								io.StringIO(mut_list[i]), 
								strict = False))
	tree_sequence_list_overall += tree_sequence_list

	print("Successfully constructed tree sequences from tables")

	# Adds neutral mutations onto the resulting tree sequence
	# This is done for each tree sequence within the list generated above, all with 
	# the same neutral mutation rate and model
	rate = 0.00001
	for i in range(len(tree_sequence_list)):
		tree_sequence_list[i] = msprime.sim_mutations(tree_sequence_list[i], rate = rate, model = msprime.InfiniteAlleles(), discrete_genome = False, keep = False)

	print("Successfully simulated neutral mutations")

	# Generates the unfolded site frequency spectrum for each tree within the list
	site_frequency_spectrum_list = []
	for i in range(len(tree_sequence_list)):
		site_frequency_spectrum_list.append(tree_sequence_list[i].allele_frequency_spectrum(span_normalise = False, polarised = True))
	print("Successfully generated site frequency spectrum")

	# Finally, we average these site frequency spectra to form a single time-averaged
	# site frequency spectrum for the run currently being considered 
	siteFrequencySpectrum = []
	for i in range(len(site_frequency_spectrum_list[0])):
		cur_sum = 0
		for j in range(len(site_frequency_spectrum_list)):
			cur_sum += site_frequency_spectrum_list[j][i]
		cur_avg = cur_sum / len(site_frequency_spectrum_list)
		siteFrequencySpectrum.append(cur_avg)
		
	# This site frequency spectrum is set to the appropriate external variable
	if link_num == 500:
		siteFrequencySpectrum_500 = siteFrequencySpectrum
	elif link_num == 1000:
		siteFrequencySpectrum_1000 = siteFrequencySpectrum
	elif link_num == 1500:
		siteFrequencySpectrum_1500 = siteFrequencySpectrum
	else:
		siteFrequencySpectrum_2000 = siteFrequencySpectrum
    

# Defines the inverse function f(x) = 1 / x for convenience in plotting
def invFunc(x):
    return 1 / x

# Defines a helper function for the Cvijovic et al. distortions
def g(x, n, ud, sd):
    return 1 / (1 + np.exp(n * sd * x * np.exp(-ud / sd)))

# Defines a piecewise function according to the distortions outlined by Cvijovic et al. (2018)
def distortions(x, n, sigma, ud, sd, un, expon):
      if (x < 1 / (N * sigma)):
            return (2 * n * un) / x
      elif (x > 1 / (n * sigma) and x < expon / (N * sd)):
            return (n * un) / (n * sd * x * x * np.sqrt((ud / sd) * np.emath.logn((ud / sd), (expon / (n * sd * x)))))
      elif (x > expon / (n * sd) and x < 1 - (expon / (n * sd))):
            return (2 * n * un) / (expon * x)
      elif (x > 1 - (expon / (n * sd)) and x < 1 - (1 / (n * sigma))):
            return (n * un) / (n * sd * (1 - x) * np.emath.logn((ud / sd), (expon / (n * sd * (1 - x)))))
      elif (x > 1 - (1 / (n * sigma))):
            return (2 * n * un) / x
      else:
            return np.nan

    
# Sets up any necessary parameters, which are determined via keyboard input for ease in updating
N = 1000
Sd = 0.009487
# To find theta, we average the values for our tree sequences
theta_sum = 0
for i in range(len(tree_sequence_list_overall)):
     theta_sum += tree_sequence_list_overall[i].diversity(mode = "branch")
theta = theta_sum / len(tree_sequence_list_overall)
coal_ne = theta / (4 * N)

# Sets the number of bars to consider for the histogram (in case only a partial SFS is desired)
bars = len(siteFrequencySpectrum_500)

# Constructs an appropriate x-axis based on Cvijovic et al.
X = [x / (2 * N) for x in range(1, bars)]

# Sets appropriate axis labels for the overall plot
plt.xlabel("Allele count / 2N", fontweight = "bold", fontsize = 20)
plt.ylabel("Normalized number of sites", fontweight = "bold", fontsize = 20)

# Constructs the arrays for both of the neutral curves, along with for the Cvijovic et al. distortions
inverses_n = []
inverses_coal_ne = []
for i in range(0, len(X)):
    inverses_n.append((2 * N * rate) / X[i])
    inverses_coal_ne.append((2 * coal_ne * rate) / X[i])
    
# Cleans the distortions data, removing any significant spikes and replacing 
# them simply with NaN values 
# Note that the adjustments to the upper and lower bound will likely need to be adjusted
# for different parameter ranges
'''
upper_bound = 1 - (pos_expon / (N * Sd)) + 0.1
lower_bound = pos_expon / (N * Sd) - 0.11
for i in range(0, len(X)):
	if (distortion_list[i] > inverses_n[i]):
		distortion_list[i] = np.nan
	elif (X[i] > lower_bound and X[i] < upper_bound and distortion_list[i] != inverses_cvijovic_ne[i]):
		distortion_list[i] = np.nan
'''

    
# Normalizing the neutral curves & Cvijovic et al. distortions
norm_val = inverses_n[0]
for i in range(0, len(X)):
    inverses_n[i] /= norm_val
    inverses_coal_ne[i] /= norm_val
    
# Normalizes the site frequency spectra to better match Cvijovic et al. figures
norm_SFS_500 = []
norm_SFS_1000 = []
norm_SFS_1500 = []
norm_SFS_2000 = []
total_sum_500 = sum(siteFrequencySpectrum_500)
total_sum_1000 = sum(siteFrequencySpectrum_1000)
total_sum_1500 = sum(siteFrequencySpectrum_1500)
total_sum_2000 = sum(siteFrequencySpectrum_2000)
for i in range(0, len(siteFrequencySpectrum_500)):
	norm_SFS_500.append(siteFrequencySpectrum_500[i] / (norm_val * 11500))
	norm_SFS_1000.append(siteFrequencySpectrum_1000[i] / (norm_val * 23000))
	norm_SFS_1500.append(siteFrequencySpectrum_1500[i] / (norm_val * 34500))
	norm_SFS_2000.append(siteFrequencySpectrum_2000[i] / (norm_val * 46000))
    
# Uses a moving average to try and "smooth" the actual SFS, with this being done for
# each of the possible quantities of linkage blocks
clean_SFS_500 = pd.DataFrame({"data" : norm_SFS_500[1:]})
clean_SFS_500 = clean_SFS_500.rolling(window = 50).mean()
clean_SFS_500 = clean_SFS_500['data'].to_list()
i = 0
while (np.isnan(clean_SFS_500[i])):
     clean_SFS_500[i] = norm_SFS_500[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_500[i + j] = norm_SFS_500[i + j]
	 
clean_SFS_1000 = pd.DataFrame({"data" : norm_SFS_1000[1:]})
clean_SFS_1000 = clean_SFS_1000.rolling(window = 50).mean()
clean_SFS_1000 = clean_SFS_1000['data'].to_list()
i = 0
while (np.isnan(clean_SFS_1000[i])):
     clean_SFS_1000[i] = norm_SFS_1000[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_1000[i + j] = norm_SFS_1000[i + j]
	 
clean_SFS_1500 = pd.DataFrame({"data" : norm_SFS_1500[1:]})
clean_SFS_1500 = clean_SFS_1500.rolling(window = 50).mean()
clean_SFS_1500 = clean_SFS_1500['data'].to_list()
i = 0
while (np.isnan(clean_SFS_1500[i])):
     clean_SFS_1500[i] = norm_SFS_1500[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_1500[i + j] = norm_SFS_1500[i + j]
	 
clean_SFS_2000 = pd.DataFrame({"data" : norm_SFS_2000[1:]})
clean_SFS_2000 = clean_SFS_2000.rolling(window = 50).mean()
clean_SFS_2000 = clean_SFS_2000['data'].to_list()
i = 0
while (np.isnan(clean_SFS_2000[i])):
     clean_SFS_2000[i] = norm_SFS_2000[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_2000[i + j] = norm_SFS_2000[i + j]
     
# Alternatively, uses a Savitzky-Golay filter in order to smooth particularly noisy data
'''
golay_SFS = savgol_filter(norm_SFS[1:], 200, 3, mode="nearest")
for i in range(0, len(golay_SFS)):
     if golay_SFS[i] <= 10 ** -9:
          golay_SFS[i] = np.nan
'''
          
# Plots the actual SFS, ensuring agreement with current axis length for each linkage block quantity
while(len(clean_SFS_500) != len(X)):
     norm_SFS_500.append(np.nan)
     clean_SFS_500.append(np.nan)
while(len(clean_SFS_1000) != len(X)):
     norm_SFS_1000.append(np.nan)
     clean_SFS_1000.append(np.nan)
while(len(clean_SFS_1500) != len(X)):
     norm_SFS_1500.append(np.nan)
     clean_SFS_1500.append(np.nan)
while(len(clean_SFS_2000) != len(X)):
     norm_SFS_2000.append(np.nan)
     clean_SFS_2000.append(np.nan)
# plt.plot(X, norm_SFS[1:], "o", label = "Observed", color = "#ffcaca")
          
# Plotting both neutral curves for N and Ne
# plt.plot(X, inverses_n, label = "Neutral with N", color = "#f47a00", linewidth = 5)
# plt.plot(X, inverses_cvijovic_ne, label = "Cvijovic Reduced Ne", color = "darkorchid")
# plt.plot(X, inverses_joseph_ne, label = "Joseph Reduced Ne", color = "deeppink") #62c8d3 
plt.plot(X, inverses_coal_ne, label = "Neutral with Ne", color = "#073b4c", linewidth = 6)
    
# Plots the cleaned site frequency spectra, separately for each possible number of linkage blocks
plt.plot(X, clean_SFS_500, label = "Moving average with L = 11,500", color = "#118ab2", linewidth = 3, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_1000, label = "Moving average with L = 23,000", color = "#06d6a0", linewidth = 4, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_1500, label = "Moving average with L = 34,500", color = "#ffd166", linewidth = 5, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_2000, label = "Moving average with L = 46,000", color = "#ef476f", linewidth = 6, linestyle = "dashed", dashes = (3, 3))
# plt.plot(X, golay_SFS, label = "Savitzky-Golay Smoothed SFS", color = "r", linewidth = 2)

# Plotting the Cvijovic et al. distortions
# plt.plot(X, distortion_list, label = "Cvijovic et al.", color = "#007191", linewidth = 5, linestyle = "dashed", dashes = (3, 3))

# Sets the scaling, legend, and ensures appropriate structure of the plot
plt.xscale("logit")
plt.yscale("log")
plt.ylim(10 ** -5, 2.0)

# Adjusts the order in the legend to improve readability
handles, labels = plt.gca().get_legend_handles_labels()
order = [1, 2, 0, 3, 4]
plt.legend([handles[i] for i in order], [labels[i] for i in order], fontsize = 16)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.tight_layout()

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(13)
plt.savefig("Block_Comparison_Plot.png")