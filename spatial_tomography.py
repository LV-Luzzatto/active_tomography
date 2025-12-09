'''
Minimal code to implement spatial cluster tomography.

A version of this code was used for parts of the analysis in the paper "Spatial and Temporal
Cluster Tomography of Active Matter", by L.V.Luzzatto, M. Casiulis, S. Martiniani, & I.A. Kovács,
available as an e-print at https://arxiv.org/abs/2511.09444.

Requires: cluster_funcions.py

Leone Luzzatto - created: 2023.11.20; updated: 2024.11.15; final version: 2025.11.25.
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
'''
def bin_positions(positions, size):
    occupied = np.zeros([size, size])
    for pos in positions:
        x = int(pos[0])
        y = int(pos[1])
        occupied[x,y] = 1
    return occupied


#__________________________________________
'''
analysis(binned, size, LC='exclude')
Run the cluster tomography analysis on a binned configuration of the particles

Input:
    - (array) [size,size] array containing the binned particle configuration from one snapshot of the data
    - (int) linear size of the system (in units of the particle diameter)

Options:
    - how to treat the largest cluster (LC): 'include' or 'exclude'

Return:
    - corner contribution, C(ell) (float)
    - gap-size statistics, g(s) (array)
    - size of the largest cluster (int)
'''
def analysis(binned, size, LC='exclude'):
    
    # identify contiguous clusters from the binned configuration w/ Hoshen-Kopelman algorithm
    bonds = cf.get_bonds_nd(binned, d=2, bond_type='one')
    labels = cf.hoshen_kopelman_nd(bonds, [size,size], d=2, bc='pp')
    
    # LC stands for "largest cluster"
    LC_label, LC_size = Counter(labels).most_common(1)[0][:]

    # exclude largest cluster (if needed)
    if LC=='exclude':
        labels = np.where( labels==LC_label, np.arange(size*size), labels )

    # gap-size statistics, g(s)
    gap_stat = cf.get_gapstat_2d(labels, size, bc='p', both_directions=True)
    
    # corner contribution, C(ell)
    corner = np.dot( np.arange(size//2), gap_stat[:size//2] + np.flip(gap_stat[size//2:]) ) / (2 * size)

    return corner, gap_stat, LC_size



####################################
# ANALYSIS

t_start = perf_counter()

# Set variables
size = 64
Pe_r = 10
phi = 0.515
snapshots = 100
LC = 'exclude'

# Set up input and output files
from pathlib import Path
in_dir = './example_spatial_data/Pe_r{}/L{}/Phi{:.3f}/'.format(Pe_r,size,phi)
out_dir = './spatial_tomography_results/Pe_r{}/L{}/Phi{:.3f}/'.format(Pe_r,size,phi)
Path(out_dir).mkdir(parents=True, exist_ok=True)

# lists to hold the results
cc = [] # corner contribution
gs = [] # gap-size statistics
lc = [] # size of the largest cluster

# Analyze the data
for dump in range(1, snapshots+1):
    with open(in_dir + 'dump.{:08d}.txt'.format(dump), 'r') as data_file:
        positions = np.loadtxt(data_file, skiprows=2, usecols=(0,1))[:,:2]
    
    binned = bin_positions(positions, size)
    corner, gap_stat, LC_size = analysis(binned, size, LC=LC)

    cc.append(corner)
    gs.append(gap_stat)
    lc.append(LC_size)

# Write results to .txt files
with open(out_dir + 'corner.txt', 'w') as outfile: 
    np.savetxt(outfile, cc, header='Corner contribution, C(ell)', delimiter='\n')
with open(out_dir + 'gapstat.txt', 'w') as outfile: 
    np.savetxt(outfile, np.array(gs), header='Gap-size statistics, g(s)', delimiter='\t')
with open(out_dir + 'LC.txt', 'w') as outfile: 
    np.savetxt(outfile, lc, header='Size of the largest cluster', delimiter='\n')


###################################
# END

t_stop = perf_counter()
print('Spatial cluster tomography')
print('L = {}; Pe_r = {}; phi = {:.3f}; {} snapshots'.format(size, Pe_r, phi, snapshots))
print('Results saved ({:.2f} seconds).'.format(t_stop-t_start))


