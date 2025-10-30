#!/usr/bin/python

'''
The following program builds on the tskit library in order to analyze the succinct
tree sequences constructed by the simulation and to generate a site frequency spectrum
from these. Note that this code requires numpy, matplotlib, and msprime to be included
in the current working environment. Note that this version of the code only generates
the SFS, along with outputting Tajima's D.

Author: Micaila Marcelle
'''

# Imports all necessary libraries/packages
import math
import tskit
import io
import numpy as np
import matplotlib.pyplot as plt
import msprime

print("Successfully imported necessary modules")

# Reads in the necessary tskit tables 
with open("nodetable.txt") as file:
    node_table = file.read()
with open("edgetable.txt") as file:
    edge_table = file.read()
with open("sitetable.txt") as file:
    site_table = file.read()
with open("mutationtable.txt") as file:
    mut_table = file.read()
    
print("Successfully opened table files")
    
# Uses these tables in order to construct the corresponding tree sequence
# Note: once this overall file is working, as Ulises mentioned, it would be
# a good idea to do some experimentation for whether the sites/mutations 
# tables really need to be included.
treeSequence = tskit.load_text(io.StringIO(node_table), 
                               io.StringIO(edge_table), 
                               io.StringIO(site_table),
                               io.StringIO(mut_table),
                               strict = False)

print("Successfully constructed tree sequence from tables")

# Adds neutral mutations onto the resulting tree sequence
# Note: this rate should be changed later! Currently just a filler value- should
# be adjusted for greater accuracy related to program needs
treeSequence = msprime.sim_mutations(treeSequence, rate = 0.000001, model = msprime.InfiniteAlleles(), discrete_genome = True, keep = False)

print("Successfully simulated neutral mutations")

# Generates the site frequency spectrum
# Note that this is currently obtaining the unfolded site frequency spectrum
siteFrequencySpectrum = treeSequence.allele_frequency_spectrum(span_normalise = False, polarised = True)
print("Successfully generated site frequency spectrum")

# Uses this site frequency spectrum to generate the associated histogram
# Likely will need to change the labels on the histogram? Not sure exactly what they should be
# Note that this code is heavily derived from http://sesame.uoregon.edu/~adkern/stdpopsim/doc/tutorial.html 
# To test that analytical expectations are matched, optional adjustments have also been made to create
# pairs of bars, with the latter expressing those expectations
bars = 21
expectations = []
stdErrors = []
stdErrorsMult = []
for i in range(0, bars):
    expectations.append(i * siteFrequencySpectrum[i])
    stdErrorsMult.append(i * math.sqrt(siteFrequencySpectrum[i]))
    stdErrors.append(math.sqrt(siteFrequencySpectrum[i]))

X = [x for x in range(0, bars)]
r = np.arange(bars)
plt.xlabel("Allele count", fontweight = "bold")
plt.ylabel("Number of sites", fontweight = "bold")
plt.bar(r - 0.15, yerr = stdErrors, height = siteFrequencySpectrum[0:bars], width = 0.3, label = "Actual SFS", color = "r", edgecolor = "black")
plt.bar(r + 0.15, yerr = stdErrorsMult, height = expectations, width = 0.3, label = "SFS multiplied by i", color = "c", edgecolor = "black")
plt.xticks(r, X)
plt.tight_layout()
plt.legend()
plt.savefig("alleleFrequencySpectrum.png")

# Finally, Tajima's D is determined and printed to the terminal
print("Tajima's D: " + str(treeSequence.Tajimas_D()))