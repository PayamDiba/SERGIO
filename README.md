# SERGIO (Single-cell ExpRession of Genes In silico)
SERGIO v1.0.0

Saurabh Sinha’s Lab, University of Illinois at Urbana-Champaign [Sinha Lab](https://www.sinhalab.net/sinha-s-home)

Developed by Payam Dibaeinia

## Description
SERGIO is a simulator for single-cell expression data guided by gene regulatory networks. A command-line, easy-to-use version of SERGIO will be soon uploaded to PyPI. Here is the documentation for using SERGIO v1.0.0 as a module in python.

## Dependencies
Python >= 2.7.14

numpy

scipy

networkx


## Usage
1. An instance of SERGIO simulator can be constructed as below:

```python
from sergio import sergio
sim = sergio(number_genes, number_bins, number_sc, noise_params,
    noise_type, decays, dynamics, sampling_state, dt, bifurcation_matrix, 
    noise_params_splice, noise_type_splice, splice_ratio, dt_splice)
```

* number_genes: total number of genes present in GRN
* number_bins: total number of distinct cell types to be simulated
* number_sc: total number of cells per cell type to be simulated
* noise_params: a single scalar or a list of size number_genes containing the genes’ noise amplitude parameter q in steady-state simulations or unspliced transcripts’ noise parameter in differentiation simulations. For differentiation simulations, small values (<0.5) are recommended.
* noise_type: The type of genes' stochastic noise in steady-state or unspliced transcripts' noise type in differentiation simulations. Options: “dpd”, “sp”, “sd” (For more details, see the paper)
* decays: a single scaler or a list of size number_genes containing the genes’ decay parameter for steady-state simulations or unspliced transcripts’ decay in differentiation simulations.
* sampling_state: an integer determining the length of simulations in stationary region. In steady-state simulations, for each cell type, simulations are continued for sampling_state times number_sc time steps after reaching to steady-state region. In differentiation simulations, for each cell type, if takes n steps till reaching to steady-state, simulations are continued for sampling_state times n more time steps in steady-state region. 
* dt: integration time step in steady-state simulations (default: 0.01).
* dynamics: a Boolean showing whether to simulate steady-state (False) or differentiation (True).
* bifurcation_matrix: only needed for dynamics simulations (default: None). A 2d (number_bins times number_bins) python list containing >=0 floats showing the differentiation graph. The element in row i and column j shows the migration rate (r) from cell type i to cell type j. Therefore, r times number_sc paths between cell type i and j is simulated. Increasing r slows down simulations but increases the density of simulated cells differentiating from cell type i to j, also r=0 denotes no differentiation from cell type i to j. Typically values of r around 1 result in desirable differentiation trajectories.
	- Example: system of three cell types with linear differentiation graph:
	    bifurcation_matrix = [[0, 0.8, 0 ],[0, 0, 1.1], [0,0,0]]

* noise_params_splice: only needed for dynamics simulations (default: None). A single scalar or a list of size number_genes containing the spliced transcripts’ noise parameter. Small values (<0.5) are recommended.
* noise_type_splice: only needed for dynamics simulations (default: None). The type of stochastic noise for simulations of spliced transcripts. Options: “dpd”, “sp”, “sd” (For more details, see the paper)
* splice_ratio: only needed for dynamics simulations (default: 4). A single scalar or a list of size number_genes containing the ratio of the expected expression of spliced to unspliced transcripts of genes in differentiation simulations. This tunes the degradation rate of spliced RNA.
* dt_splice: only needed for dynamics simulations (default: 0.01). Integration time step in differentiation simulations.
