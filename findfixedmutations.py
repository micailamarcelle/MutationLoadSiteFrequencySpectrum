import os
import io
import numpy as np
import pandas as pd
import sys
import concurrent.futures
import msprime
import tskit
import math
import matplotlib.pyplot as plt

summarydatafilename = "treesequencedataforUb_001_sb_002.txt"
summarydatafile = open(summarydatafilename, "w")


#import necessary packages
import tskit
import math
import sys
import io
import numpy as np
import pandas as pd
import os

#read tables for tskit
with open('sitetable.txt') as f:
    sites = f.read()

with open('nodetable.txt') as f:
    nodes = f.read()

with open('mutationtable.txt') as f:
    mutations = f.read()

with open('edgetable.txt') as f:
    edges = f.read()

#load in the tree sequence data
ts = tskit.load_text(
    nodes = io.StringIO(nodes),
    edges = io.StringIO(edges),
    sites = io.StringIO(sites),
    mutations = io.StringIO(mutations),
    strict = False)

#ts_2 = TableCollection.tree_sequence("tables.trees")

num_samples = ts.get_sample_size()
print("the size of the sample of text-based ts is " + str(num_samples) + ".")
#num_samples = ts_2.get_sample_size()
#print(f"the size of the sample of direct load ts is {num_samples}.")

N = int(input("Enter popsize:"))
burnin = float(input("Enter burn-in generation:"))
G = float(N*10)

threshold_gen = (G-burnin)*N
#output the fixation time table for fixed mutations in final sample
with open('fixed_mut.txt', 'w', encoding='utf-8') as f:
    print('mutation', sep='\t', file=f)
    for tree in ts.trees(root_threshold=2):
        for root in tree.roots:
            # Create a list to store mutations that meet the time threshold
            filtered_mutations = []
            for mutation in tree.mutations():
                # Get the associated node ID for the mutation
                node_id = mutation.node

                # Access the time attribute of the node
                node_time = tree.time(node_id)
                
                if node_time <= threshold_gen:
                    filtered_mutations.append(mutation)
            
            # Iterate through the filtered mutations to find those where root == mutation.node
            for mutation in filtered_mutations:
                if root == mutation.node:
                    print(mutation.derived_state, sep='\t', file=f)

mutfixtable = pd.read_csv('fixed_mut.txt', sep="\t", header=0)

mut_ben = mutfixtable[mutfixtable["mutation"] > 0]
mut_del = mutfixtable[mutfixtable["mutation"] < 0]
#printing
print("Beneficial flux is: ", f'{((2 * (mut_ben["mutation"] + 1).apply(math.log)).sum())/(G-burnin):.10f}')
print("Deleterious flux is: ", f'{((2 * (mut_del["mutation"] + 1).apply(math.log)).sum())/(G-burnin):.10f}')

