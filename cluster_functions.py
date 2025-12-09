'''
Functions used to identify and study the cluster structure of n-dimensional data.

A version of this code was used for parts of the analysis in the paper "Spatial and Temporal
Cluster Tomography of Active Matter", by L.V.Luzzatto, M. Casiulis, S. Martiniani, & I.A. Kovács,
available as an e-print at https://arxiv.org/abs/2511.09444.

Leone Luzzatto - created: 2024.02.21; updated: 2025.01.14; final version: 2025.11.25.
'''

import numpy as np
from itertools import product


#______________________________
'''
find(x, labels)
Required for the Hoshen-Kopelman cluster-finding algorithm.

Input:
    - (int) position in the array of labels
    - (array) cluster labels

Return:
    - new label for position x
'''
def find(x, labels):
    x = int(x)
    while labels[x] != x:
        x = labels[x]
    return x


#______________________________
'''
hoshen_kopelman_nd(bonds, sizes, d=2, bc='pp')
Hoshen-Kopelman algorithm to identify contiguous clusters in any number of dimensions.

NOTE    Any list of properties of the system along different axes are ordered as (... , z, y, x)

Input:
    - (list of arrays) connectivity matrices, one for each direction of the bonds (..., z, y, x).
        Open and closed bonds are marked as 1 and 0.
    - (list of integers) system sizes in each direction

Options:
    - (int) number of dimensions
    - boundary conditions (bc) along each direction. 'o' for Open, 'p' for Periodic.
        e.g., bc='ppo' means a 3d system with periodic bc along the z and y directions and open bc along x.
    - reduce the label list to use the smallest possible integer sequence
        e.g., [0 1 1 3 1 5 5 1] -> [0 1 1 2 1 3 3 1])
    
Return:
    - array of cluster labels
'''
def hoshen_kopelman_nd(bonds, sizes, d=2, bc='pp', reduce=False):
    labels = np.arange(np.prod(sizes)) # list containing each site's label
    convert_vec = [np.prod(sizes[k+1:]).astype(int) for k in range(d)] # helpful for converting coordinates: (... z, y, x) -> (x + size*y + size**2*z + ...)

    # coordinates in the form (... , z, y, x)
    for coords in product(*[range(size) for size in sizes]):
        for n in range(d):
            if coords[n] != 0:
                if bonds[n][coords]:
                    idx_site = np.dot(coords, convert_vec)
                    idx_neighbor = idx_site - convert_vec[n]
                    labels[idx_site] = find(idx_neighbor, labels)
                    
                    for m in range(n+1,d):
                        if coords[m] != 0:
                            if bonds[m][coords]: # activated bond in direction m and n
                                idx_neighbor = idx_site - convert_vec[m]
                                labels[find(idx_neighbor, labels)] = find(idx_site, labels)
                    
                    break
        
    # apply periodic boundary conditions
    for n in range(d):
        if bc[n] == 'p':
            positions = [range(size) for size in sizes[:n]] + [[0]] + [range(size) for size in sizes[n+1:]]
            
            for coords in product(*positions):
                if bonds[n][coords]: # activated bond in direction n
                    idx_site = np.dot(coords, convert_vec)
                    idx_neighbor = idx_site + convert_vec[n]*(sizes[n]-1)
                    
                    labels[find(idx_neighbor, labels)] = labels[find(idx_site, labels)]

    # label each size with its corresponding equivalence class (i.e. cluster)
    for idx in range(len(labels)):
        labels[idx] = int(find(idx, labels))
    
    if reduce:
        labels = np.unique(labels, return_inverse=True)[1]
    
    return labels


#______________________________
'''
get_bonds_nd(state, d=2, bc='pp', bond_type='one')
Given data on a lattice (square, cubic, hypercubic), identify the connectivity matrices.

Input:
    - (array) state of the system. The value at each site marks its state.

Options:
    - (int) number of dimensions
    - bond type, 'one' or 'all'.
        If 'one', only sites in state '1' form clusters. Use for e.g. site percolation.
        If 'all', each site connects with neighboring sites in the same state. Use for e.g. Ising model.
      
Return:
    - list of connectivity matrices, one for each direction of the bonds (..., z, y, x).
'''
def get_bonds_nd(state, d=2, bond_type='one'):
        
    rolls = [np.roll(state, 1, axis=n) for n in range(d)]
    
    if bond_type == 'all':
        bonds = [roll==state for roll in rolls]
        
    elif bond_type == 'one':
        bonds = [roll*state for roll in rolls]
    
    return bonds


#______________________________
'''
get_gapstat_nd(labels, size, bc='p', both_directions=False, bad=np.nan)
Calculate the gap-size statistics, g(s), from the list of cluster labels.

NOTE    *** This function is only for 2d data ***

Input:
    - (array) cluster labels
    - (int) linear size of the system

Options:
    - boundary conditions (bc). 'o' for Open, 'p' for Periodic.
    - whether gaps should be measured using both horizontal and vertical lines
    - bad data to be excluded
      
Return:
    - gap-size statistics, g(s)
'''
def get_gapstat_2d(labels, size, bc='p', both_directions=False, bad=np.nan):
    f = np.zeros([size,])
    norm = size
    
    completed = False
    while not completed:

        for row in range(size):
            line_labels = labels[row*size : (row+1)*size]

            first_seen = {}
            last_seen = {}

            for i, lab in enumerate(line_labels):
                
                if lab != bad:
                    
                    if lab in last_seen:
                        gap = i - last_seen[lab]
                        f[gap] += 1

                    else:
                        first_seen[lab] = i

                    last_seen[lab] = i

            if bc == 'p':
                for lab in last_seen:
                    gap = size + first_seen[lab] - last_seen[lab]
                    f[gap%size] += 1
        
        if both_directions:
            labels = np.transpose(labels.reshape([size,size])).ravel()
            norm *= 2
            both_directions = False
        else:
            completed = True
    
    return f / norm


