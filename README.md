# Active Tomography
Minimal Python implementation of spatial and temporal *cluster tomography* on snapshots of a two-dimensional active system.

A version of this code was used for parts of the analysis in the paper "Spatial and Temporal
Cluster Tomography of Active Matter," by L.V. Luzzatto, M. Casiulis, S. Martiniani, & I.A. Kovács,
available as an e-print at https://arxiv.org/abs/2511.09444, where this analysis is applied to simulations of active Brownian Particles (see https://github.com/mkasiulis/active_brawnian for the simulation code).


## Description
The code is split between four files:
1. `cluster_functions.py` contains several functions needed to identify different clusters and measure relevant quantities.
2. `spatial_tomography.py` implements spatial cluster tomography on a collection of snapshots.
3. `cluster_membership.py` identifies clusters from a collection of snapshots and saves the results in a format that can be used by `temporal_tomography.py`.
4. `temporal_tomography.py` implements temporal cluster tomography on the output of `cluster_membership.py`.

Two small example datasets are included to test the code.

### `cluster_functions.py`
Most functions contained in this file have optional arguments so they can be used for different systems, including:
 - $n$-dimensional systems, with $n\geq2$,
 - periodic, open, or any combination of periodic and open boundary conditions.

### `spatial_tomography.py`
implement spatial cluster tomography in a square system with periodic boundary conditions. 

A few variables need to be specified in the code:
 - the linear size of the system,
 - the number of snapshots to be analyzed.

Additionally, the volume fraction, $\phi$, and angular Péclet number, $\mathrm{Pe}_r$, need are required to locate the example data ($\phi=0.515$ and $\mathrm{Pe}_r=10$).

The code will create three separate output files containing the results of the analysis (see [the paper](https://arxiv.org/abs/2511.09444) for the relevant definitions):
 1. `corner.txt` contains the *corner contribution*, $\mathcal{C}$, calculated from each sample,
 2. `gapstat.txt` contains the *gap-size statistics*, $g(s)$,
 3. `LC.txt` contains the volume of the largest cluster in each sample.

 > **NOTE.** Boundary conditions should be handled carefully. If a system has periodic boundaries, a naive implementation typically returns a large number of small gaps, as well as a large number of large gaps, with size close to the linear size of the system, due to the periodic boundaries. In this case, it is often useful to account for the periodic boundaries by redefining gap sizes as $\tilde{s}=\mathrm{min}$($s,L-s$), with $L$ the linear size of the system, and study $g$($\tilde{s}$). See, for example, appendix B and the SM of "Cluster Tomography in Percolation," H. Ansell, S. Frank, & I.A. Kovács, Phys. Rev. Research 5, 043218 (2023), https://arxiv.org/abs/2307.04260.



### `cluster_membership.py`
Identify the cluster membership of each particle in the system.

Variables to be specified in the code:
 - linear size of the system,
 - number of snapshots to be analyzed,
 - volume fraction, $\phi=0.26$, and angular Péclet number, $\mathrm{Pe}_r=100$ (needed to locate the example data).

The code creates two output files:
 1. `clusters.out` is a binary file containing the cluster membership of each particle in the system in each snapshot. See the code for a description of the format.
 2. `LC.txt` contains the volume of the largest cluster in each snapshot.


### `temporal_tomography.py`
Implement temporal cluster tomography from cluster membership data. Designed to run on the output of `cluster_membership.py`. 

Variables to be specified in the code:
 - linear size of the system,
 - number of time-steps (i.e., snapshots) to be analyzed,
 - number of particles used in the analysis,
 - volume fraction, $\phi=0.26$, and angular Péclet number, $\mathrm{Pe}_r=100$ (needed to locate the example data).

 The analysis considers all particle pairs within a subset of the total number of particles. Each particle $i$ is assigned a unique numerical ID across all snapshots, $n_i = 0,1,2,3,...$ . For ease of implementation, given a number $N$ of particles to be used in the analysis, the code considers all particles whose IDs are $n_i<N$. Since (1) the particles are initialized in random positions and (2) they mix as the system reaches a steady-state, this selection is effectively a random sample of $N$ particles.

 The code creates the output file `temporal_gapsize.txt`, that contains the average temporal gap-size statistics $g_t(s_t)$ for the specified particles across the selected snapshots.

 > **NOTE.** Numba is used to speed up the code, which requires some data types to be specified for the purpose of memory allocation. The chosen data types (e.g. int64, uint32, bool, etc.) might have to be modified to suit larger/smaller numbers of particles.
