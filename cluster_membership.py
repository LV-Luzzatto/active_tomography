'''
Minimal code to identify cluster membership of particles over time. Needed to perform temporal
cluster tomography.

A version of this code was used for parts of the analysis in the paper "Spatial and Temporal
Cluster Tomography of Active Matter", by L.V.Luzzatto, M. Casiulis, S. Martiniani, & I.A. Kovács,
available as an e-print at https://arxiv.org/abs/2511.09444.

Requires: cluster_funcions.py

    
NOTE    One of the output files, 'clusters.out' is a binary pickle file containing the cluster
        membership of different particles across snapshot. It can be read using 'pickle.load()'.
        It is formatted as a list of particle IDs in which '-1' is used to separate
        different clusters.
        For example, the list [0, 1, 4, -1, 2, 5, -1, 7, 8, 10, 11] means that at that snapshot
        there are three clusters: [0, 1, 4], [2, 5], and [7, 8, 10, 11], where each non-negative
        number uniquely identifies one particle across snapshots.


Leone Luzzatto - created: 2025.01.14; final version: 2025.11.25.
'''

import numpy as np
import cluster_functions as cf
from time import perf_counter
from collections import Counter


###########################################
# FUNCTIONS

#__________________________________________
'''
bin_positions(positions, size)
Bin the particle positions.

Input:
    - (array) position (x,y) of the center of each particle (in units of the particle diameter)
    - (int) linear size of the system (in units of the particle diameter)

Return:
    - array of occupied and unoccupied bins, marked by 1 or 0
    - array with the index of the bin to which each particle belongs
'''
def bin_positions(positions, size):
    occupied = np.zeros([size, size])
    particle_bin = np.zeros([len(positions),], dtype=np.int_)
    for i, pos in enumerate(positions):
        x = int(pos[0])
        y = int(pos[1])
        occupied[x,y] = 1
        particle_bin[i] = size*y + x
    
    return occupied, particle_bin

#_________________________________________________________________________________________________________
'''
cluster_membership(binned, size)
Identify the cluster membership of each particle in the system.

Input:
    - (array) [size,size] array containing the binned particle configuration from one snapshot of the data
    - (array) index of the bin to which each particle belongs (bin (x,y) -> index x+(y*size))
    - (int) linear size of the system (in units of the particle diameter)

Return:
    - dictionary containing the list of particles belonging to each cluster
    - size of the largest cluster (int)
'''
def cluster_membership(binned, particle_bin, size):
    bonds = cf.get_bonds_nd(binned, d=2, bond_type='one')
    labels = cf.hoshen_kopelman_nd(bonds, [size,size], d=2, bc='pp')
    LC_size = Counter(labels).most_common(1)[0][1]

    # make a dictionary containing the list of particles belinging to each cluster
    particle_labels = labels[particle_bin]
    clusters = [[] for _ in range(np.max(labels)+1)]
    for i, lab in enumerate(particle_labels):
        clusters[lab].append(i)
    
    # remove all clusters with only one particle. They don't contribute to the analysis.
    to_pop = []
    for n in range(np.max(labels)):
        if len(clusters[n]) <= 1:
            to_pop.append(n)
    for n in reversed(to_pop):
        clusters.pop(n)

    return clusters, LC_size


###################################
# ANALYSIS

t_start = perf_counter()

# Set variables
size = 64
Pe_r = 100
phi = 0.26
snapshots = 100

# Set up input and output files
from pathlib import Path
in_dir = './example_temporal_data/'
out_dir = './'
Path(out_dir).mkdir(parents=True, exist_ok=True)

# lists to hold the results
cm = [] # cluster membership of each particle
lc = [] # size of the largest cluster

# Analyze the data
for dump in range(1000, 1000+snapshots):
    with open(in_dir + 'dump.{:08d}.txt'.format(dump), 'r') as data_file:
        positions = np.loadtxt(data_file, skiprows=2, usecols=(0,1))[:,:2]
    
    binned, particle_bin = bin_positions(positions, size)
    clusters, LC_size = cluster_membership(binned, particle_bin, size)
    lc.append(LC_size)

    # record the cluster membership of each particle
    cm.append([])
    for cluster in clusters[:-1]:
        for i in cluster:
            cm[-1].append(i)
        cm[-1].append(-1)
    for i in clusters[-1]:
        cm[-1].append(i)

# Write results to .txt files
import pickle
with open(out_dir + 'clusters.out', 'wb') as outfile: 
    pickle.dump(cm, outfile)
with open(out_dir + 'LC.txt', 'w') as outfile: 
    np.savetxt(outfile, lc, header='Size of the largest cluster', delimiter='\n')

#################################
# END

t_stop = perf_counter()
print('L = {}; Pe_r = {}; phi = {:.2f}; {} snapshots'.format(size, Pe_r, phi, snapshots))
print('Results saved ({:.2f} seconds).'.format(t_stop-t_start))


