#!/usr/bin/python

'''
This program is designed to construct a plot comparing the site frequency spectra
between a "base case" and a five-fold increase in Ud and Sd, both separately and 
simultaneously. With such a plot, we can better understand how Ud and Sd interact
to distort the site frequency spectrum in the case of sexual populations. 

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

print("Successfully imported necessary modules")

list_Uds = [2, 10]
list_sds = ["0.004", "0.020"]

siteFrequencySpectrum_2_004 = []
siteFrequencySpectrum_2_02 = []
siteFrequencySpectrum_10_004 = []
siteFrequencySpectrum_10_02 = []

# Used for finding an "averaged out" Ne, since this will allow for a single baseline
# to compare the distortions to
tree_sequence_list_overall = []

for cur_Ud in list_Uds:
      
	for cur_sd in list_sds: 
		# Reads in the necessary tskit tables for all generations of interest, adding
		# the results to the corresponding list for each type of table
		node_list = []
		edge_list = []
		site_list = []
		mut_list = []

		# First set reads in the tskit tables for generation 500000
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/nodetable_gen500000.txt") as file:
			node_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/edgetable_gen500000.txt") as file:
			edge_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/sitetable_gen500000.txt") as file:
			site_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/mutationtable_gen500000.txt") as file:
			mut_list.append(file.read())
			
		# Second reads in the tskit tables for generation 1000000
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/nodetable_gen1000000.txt") as file:
			node_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/edgetable_gen1000000.txt") as file:
			edge_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/sitetable_gen1000000.txt") as file:
			site_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/mutationtable_gen1000000.txt") as file:
			mut_list.append(file.read())
			
		# Third reads in tskit tables for generation 1500000
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/nodetable_gen1500000.txt") as file:
			node_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/edgetable_gen1500000.txt") as file:
			edge_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/sitetable_gen1500000.txt") as file:
			site_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/mutationtable_gen1500000.txt") as file:
			mut_list.append(file.read())
			
		# Finally, we read in tskit tables for generation 2000000
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/nodetable.txt") as file:
			node_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/edgetable.txt") as file:
			edge_list.append(file.read())  
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/sitetable.txt") as file:
			site_list.append(file.read())
		with open("datafor_relative_tskitstatus_ON_BURNIN_fixationcalc_OFF_Sb_" + cur_sd + "0" + "_deldist_point_bendist_exponential_mub_0.0000_chromnum_23_N0_1000_mud_" + str(cur_Ud) + "_L_1000_seed_24_Sd_" + cur_sd + "000" + "/mutationtable.txt") as file:
			mut_list.append(file.read())
			
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
		if  cur_Ud == 2 and cur_sd == "0.004":
			siteFrequencySpectrum_2_004 = siteFrequencySpectrum
		elif cur_Ud == 2 and cur_sd == "0.020":
			siteFrequencySpectrum_2_02 = siteFrequencySpectrum
		elif cur_Ud == 10 and cur_sd == "0.004":
			siteFrequencySpectrum_10_004 = siteFrequencySpectrum
		else:
			siteFrequencySpectrum_10_02 = siteFrequencySpectrum
    

# Defines the inverse function f(x) = 1 / x for convenience in plotting
def invFunc(x):
    return 1 / x
    
# Sets up any necessary parameters, which are determined via keyboard input for ease in updating
N = 1000
# To find theta, we average the values for our tree sequences
theta_sum = 0
for i in range(len(tree_sequence_list_overall)):
     theta_sum += tree_sequence_list_overall[i].diversity(mode = "branch")
theta = theta_sum / len(tree_sequence_list_overall)
coal_ne = theta / (4 * N)

# Sets the number of bars to consider for the histogram (in case only a partial SFS is desired)
bars = len(siteFrequencySpectrum_2_004)

# Constructs an appropriate x-axis based on Cvijovic et al.
X = [x / (2 * N) for x in range(1, bars)]

# Sets appropriate axis labels for the overall plot
plt.xlabel("Allele count / 2N", fontweight = "bold", fontsize = 20)
plt.ylabel("Normalized number of sites", fontweight = "bold", fontsize = 20)

# Constructs the arrays for both of the neutral curves
inverses_n = []
inverses_coal_ne = []
for i in range(0, len(X)):
    inverses_n.append((2 * N * rate) / X[i])
    inverses_coal_ne.append((2 * coal_ne * rate) / X[i])

    
# Normalizing the neutral curves
norm_val = inverses_n[0]
for i in range(0, len(X)):
    inverses_coal_ne[i] /= norm_val
    
# Normalizes the site frequency spectra to better match Cvijovic et al. figures
L = 23000  # Using the same number of linkage blocks in all cases
norm_SFS_2_004 = []
norm_SFS_2_02 = []
norm_SFS_10_004 = []
norm_SFS_10_02 = []
total_sum_2_004 = sum(siteFrequencySpectrum_2_004)
total_sum_2_02 = sum(siteFrequencySpectrum_2_02)
total_sum_10_004 = sum(siteFrequencySpectrum_10_004)
total_sum_10_02 = sum(siteFrequencySpectrum_10_02)
for i in range(0, len(siteFrequencySpectrum_2_004)):
	norm_SFS_2_004.append(siteFrequencySpectrum_2_004[i] / (norm_val * L))
	norm_SFS_2_02.append(siteFrequencySpectrum_2_02[i] / (norm_val * L))
	norm_SFS_10_004.append(siteFrequencySpectrum_10_004[i] / (norm_val * L))
	norm_SFS_10_02.append(siteFrequencySpectrum_10_02[i] / (norm_val * L))
    
# Uses a moving average to try and "smooth" the actual SFS, with this being done for
# each of the possible quantities of linkage blocks
clean_SFS_2_004 = pd.DataFrame({"data" : norm_SFS_2_004[1:]})
clean_SFS_2_004 = clean_SFS_2_004.rolling(window = 50).mean()
clean_SFS_2_004 = clean_SFS_2_004['data'].to_list()
i = 0
while (np.isnan(clean_SFS_2_004[i])):
     clean_SFS_2_004[i] = norm_SFS_2_004[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_2_004[i + j] = norm_SFS_2_004[i + j]
	 
clean_SFS_2_02 = pd.DataFrame({"data" : norm_SFS_2_02[1:]})
clean_SFS_2_02 = clean_SFS_2_02.rolling(window = 50).mean()
clean_SFS_2_02 = clean_SFS_2_02['data'].to_list()
i = 0
while (np.isnan(clean_SFS_2_02[i])):
     clean_SFS_2_02[i] = norm_SFS_2_02[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_2_02[i + j] = norm_SFS_2_02[i + j]
	 
clean_SFS_10_004 = pd.DataFrame({"data" : norm_SFS_10_004[1:]})
clean_SFS_10_004 = clean_SFS_10_004.rolling(window = 50).mean()
clean_SFS_10_004 = clean_SFS_10_004['data'].to_list()
i = 0
while (np.isnan(clean_SFS_10_004[i])):
     clean_SFS_10_004[i] = norm_SFS_10_004[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_10_004[i + j] = norm_SFS_10_004[i + j]
	 
clean_SFS_10_02 = pd.DataFrame({"data" : norm_SFS_10_02[1:]})
clean_SFS_10_02 = clean_SFS_10_02.rolling(window = 50).mean()
clean_SFS_10_02 = clean_SFS_10_02['data'].to_list()
i = 0
while (np.isnan(clean_SFS_10_02[i])):
     clean_SFS_10_02[i] = norm_SFS_10_02[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS_10_02[i + j] = norm_SFS_10_02[i + j]
     
          
# Plots the actual SFS, ensuring agreement with current axis length for each linkage block quantity
while(len(clean_SFS_2_004) != len(X)):
     norm_SFS_2_004.append(np.nan)
     clean_SFS_2_004.append(np.nan)
while(len(clean_SFS_2_02) != len(X)):
     norm_SFS_2_02.append(np.nan)
     clean_SFS_2_02.append(np.nan)
while(len(clean_SFS_2_02) != len(X)):
     norm_SFS_2_02.append(np.nan)
     clean_SFS_2_02.append(np.nan)
while(len(clean_SFS_2_02) != len(X)):
     norm_SFS_2_02.append(np.nan)
     clean_SFS_2_02.append(np.nan)
          
# Plotting both neutral curves for Ne
# plt.plot(X, inverses_coal_ne, label = "Neutral with Ne", color = "#073b4c", linewidth = 6)
    
# Plots the cleaned site frequency spectra, separately for each possible number of linkage blocks
plt.plot(X, clean_SFS_2_004, label = "Ud = 2, Sd = 0.004", color = "#118ab2", linewidth = 3, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_2_02, label = "Ud = 2, Sd = 0.02", color = "#06d6a0", linewidth = 4, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_10_004, label = "Ud = 10, Sd = 0.004", color = "#ffd166", linewidth = 5, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, clean_SFS_10_02, label = "Ud = 10, Sd = 0.02", color = "#ef476f", linewidth = 6, linestyle = "dashed", dashes = (3, 3))

# Sets the scaling, legend, and ensures appropriate structure of the plot
plt.xscale("logit")
plt.yscale("log")
plt.ylim(10 ** -5, 2.0)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.tight_layout()

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(13)
plt.savefig("Ud_vs_Sd_Comparison_Plot.png")