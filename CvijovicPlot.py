import math
import io
import numpy as np
import matplotlib.pyplot as plt


# Defines the functions necessary for plotting the Cvijovic et al. curves
def g(x, n, ud, sd):
    return 1 / (1 + np.exp(n * sd * x * np.exp(-ud / sd)))

def distortions(x, n, sigma, ud, sd, un, expon):
      if (x < 1 / (n * sigma)):
            return (2 * n * un) / x
      elif (x > 1 / (n * sigma) and x < expon / (n * sd)):
            return (n * un) / (n * sd * x * x * np.sqrt((ud / sd) * np.emath.logn((ud / sd), (expon / (n * sd * x)))))
      elif (x < 1 - expon / (n * sd) and x > (expon / (n * sd))):
            return (2 * n * un) / (expon * x)
      elif (x > 1 - (expon / (n * sd)) and x < 1 - (1 / (n * sigma))):
            return (n * un) / (n * sd * (1 - x) * np.emath.logn((ud / sd), (expon / (n * sd * (1 - x)))))
      elif (x > 1 - (1 / (n * sigma))):
            return (2 * n * un) / x
      else:
            return np.nan

# First, we begin by outlining the parameter values for the graph being plotted
N = int(input("Enter N: "))
Ud = float(input("Enter Ud: "))
Un = float(input("Enter Un: "))
Sd = float(input("Enter Sd: "))
sigma = np.sqrt(Ud * Sd)
pos_expon = np.exp(Ud / Sd)
Ne = N * np.exp(-Ud / Sd)

# We then enumerate the values on the x-axis within the plot
x = np.arange(10.0 ** (-5), 1.0, 0.000001)

# The corresponding y values are then generated, along with those for the neutral
# expectations for both N and Ne
cvijovic_vals = []
neutral_n = []
neutral_ne = []
for i in range(0, len(x)):
	cvijovic_vals.append(distortions(x[i], N, sigma, Ud, Sd, Un, pos_expon))
	neutral_n.append((2 * N * Un) / x[i])
	neutral_ne.append((2 * Ne * Un) / x[i])
      
# Cleans the distortions data, removing any significant spikes and replacing 
# them simply with NaN values 
# Note that the adjustments to the upper and lower bound will likely need to be adjusted
# for different parameter ranges
'''
upper_bound = 1 - (pos_expon / (N * Sd)) + 0.1
lower_bound = pos_expon / (N * Sd) - 0.11
for i in range(0, len(x)):
	if (cvijovic_vals[i] > neutral_n[i]):
		cvijovic_vals[i] = np.nan
	elif (x[i] > lower_bound and x[i] < upper_bound and cvijovic_vals[i] != neutral_ne[i]):
		cvijovic_vals[i] = np.nan
'''

     
# Normalizes the y-values within each array in such a way that maintains functional relationships
norm_val = cvijovic_vals[0]
for i in range(0, len(x)):
	if (cvijovic_vals[i] != np.nan):
		cvijovic_vals[i] /= norm_val
	neutral_n[i] /= norm_val
	neutral_ne[i] /= norm_val
     
# Finally, the results are plotted
plt.title("Cvijovic et al. Expected Distortions")
plt.xlabel("Variant Frequency (f)")
plt.ylabel("Relative Fraction of SNPs")
plt.plot(x, cvijovic_vals, label = "Cvijovic et al. Distortions", color = "r", linewidth = 3)
plt.plot(x, neutral_n, label = "Neutral N Expectation", color = "b")
plt.plot(x, neutral_ne, label = "Reduced Ne Expectation", color = "g")
plt.legend()


plt.xscale("logit")
plt.yscale("log")
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(13)
plt.savefig("Cvijovic_plot.png")