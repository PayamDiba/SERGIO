# SERGIO (Single-cell ExpRession of Genes In silicO)
[![DOI](https://zenodo.org/badge/197455649.svg)](https://zenodo.org/badge/latestdoi/197455649)

SERGIO v1.0.0


Saurabh Sinha’s Lab, University of Illinois at Urbana-Champaign [Sinha Lab](https://www.sinhalab.net/sinha-s-home)

Developed by Payam Dibaeinia

## Description
SERGIO is a simulator for single-cell expression data guided by gene regulatory networks. A command-line, easy-to-use version of SERGIO will be soon uploaded to PyPI. Here is the documentation for using SERGIO v1.0.0 as a module in python.

## Dependencies
Python >= 2.7.14

numpy >= 1.13.3

scipy >= 1.1.0

networkx >= 2.0

The tool has been succefully tested on MacOS Sierra (v10.12.6) and ScientificLinux 6.9.

## Getting Started
To download SERGIO, clone the repository via the following command (should take < 1 minute):

```git clone https://github.com/PayamDiba/SERGIO```


## Usage

run_sergio.ipynb is a jupyter notebook that runs SERGIO for steady-state and differentiation simulations as well as adding technical noise. SERGIO with an easier interface for simulations and adding technical noise will be soon uploaded to PyPI. 
### Simulating Clean Data
A synthetic data set can be simulated in four lines of python code:

1. An instance of SERGIO simulator is constructed as below:

```python
import numpy as np
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
	- Example: system of three cell types with a linear differentiation graph:
	    bifurcation_matrix = [[0, 0.8, 0 ],[0, 0, 1.1], [0,0,0]]

* noise_params_splice: only needed for dynamics simulations (default: None). A single scalar or a list of size number_genes containing the spliced transcripts’ noise parameter. Small values (<0.5) are recommended.
* noise_type_splice: only needed for dynamics simulations (default: None). The type of stochastic noise for simulations of spliced transcripts. Options: “dpd”, “sp”, “sd” (For more details, see the paper)
* splice_ratio: only needed for dynamics simulations (default: 4). A single scalar or a list of size number_genes containing the ratio of the expected expression of spliced to unspliced transcripts of genes in differentiation simulations. This tunes the degradation rate of spliced RNA.
* dt_splice: only needed for dynamics simulations (default: 0.01). Integration time step in differentiation simulations.

2. GRN structure and master regulators’ profile is fed into the simulator by invoking `build_graph` method:

```python
sim.build_graph(input_file_taregts, input_file_regs, shared_coop_state)
```
	Note: Before preparing the input files, use zero-based numerical indexing for naming all gene IDs (both master regulators and non-master regulators) in the GRN. For example if there are 10 genes in the GRN, naming them starting 0 to 9.

* input_file_taregts: path to a comma separated file containing GRN structure and its parameters. Each row in this file corresponds to a target gene in the GRN. Every row contains the parameters of the hill functions of all the regulators of that row’s target gene. 
Column order is: target gene id, number of target’s regulators, regulator ID_1,…, regulator ID_n, K_1,…,K_n, hill_coeff_1, …, hill_coeff_n

	where “K” denotes the maximum interaction strength (see equation 6 in the manuscript). For activating interactions use positive “K” and for repressive ones use negative values. Since master regulators do not have any regulator they should not be included in this file as a target gene. 
	- Example: input_file_taregets for GRN of three genes  g0 --> g1 --| g2   
	1, 1, 0, 2.5, 2   
	2, 1, 1, -1.3, 2


* input_file_regs: path to a comma separated file containing master regulators’ basal production rate in all cell types. So, if there are three cell types to be simulated, each row in this file has four entries: master regulator id, production rate cell type_1,…,  production rate cell type_3.
	- Example: input_file_regs, for GRN g0 --> g1 --| g2,  in three cell types:   
	   0, 0.5, 1.5, 3

* shared_coop_state: in case of using >0 values, the same value is used for all hill coefficients in simulations and therefore there is no need to specify these values (hill_coeff) in the input_file_taregets (they are ignored otherwise). In case of using any <=0 value, hill coefficients will be read from input_file_taregets. Recommended values of hill coefficient is between 1 and 3 (default: 0).

3. For running steady-state simulations invoke `simulate` method:
```python
sim.simulate()
```

For running differentiation simulations invole `simulate_dynamics` method:
```python
sim.simulate_dynamics()
```
4. To get the clean simulated expression matrix after steady_state simulations invoke `getExpressions` method:
```python
expr = sim.getExpressions()
```

This returns a 3d numpy array (#cell_types * #genes * #cells_per_type). To convert into a 2d matrix of size (#genes * #cells) do:
```python
expr = np.concatenate(expr, axis = 1)
```

Now each row represents a gene and each column represents a simulated single-cell. Gene IDs match their row in this expression matrix, also cell types are groupd by columns such that the first #cells_per_type columns correspond to the first simulated cell type, the next #cells_per_type columns correpond to the second cell type and ... .

To get the clean simulated expression matrix after differentiation simulations invoke `getExpressions_dynamics` method:
```python
exprU, exprS = sim.getExpressions_dynamics()
```

This returns two 3d numpy array (#cell_types * #genes * #cells_per_type) for unspliced (exprU) and spliced (exprS) transcripts. To convert them into a 2d matrix of size (#genes * #cells) do:
```python
exprU = np.concatenate(exprU, axis = 1)
exprS = np.concatenate(exprS, axis = 1)
```

Now each row represents a gene and each column represents a simulated single-cell. Gene IDs match their row in this expression matrix, also cell types are groupd by columns such that the first #cells_per_type columns correspond to the first simulated cell type, the next #cells_per_type columns correpond to the second cell type and ... .


### Adding Technical Noise
SERGIO can add three type of technical noise (outlier genes, library size, and dropouts) to the clean simulated data. These noise modules can be invoked in any combination and order. Also, there is a fourth module that converts an expression matrix to an mRNA count matrix. All of these modules work on the 3d expression matrix (not the 2d concatenated version).

First use SERGIO to simulate a clean data set and obtain the 3d expression matrix:  
In steady-state simulations:  
```python
expr = sim.getExpressions()
```

In differentiation simulations:
```python
exprU, exprS = sim.getExpressions_dynamics()
```


Here we show how to add outlier genes followed by library size and then dropouts. Please refer to the manuscript for the definitions of the input parameters to the each of the noise modules:
1. **Outlier Genes**: 

In steady-state simulations invoke the `outlier_effect` method:
```python
expr_O = sim.outlier_effect(expr, outlier_prob, mean, scale)
```

In differentiation simulations invoke the `outlier_effect_dynamics` method:
```python
exprU_O, exprS_O = sim.outlier_effect_dynamics(exprU, exprS, outlier_prob, mean, scale)
```

2. **Library Size**: 

In steady-state simulations invoke the `lib_size_effect` method:
```python
expr_O_L = sim.lib_size_effect(expr_O, mean, scale)
```

In differentiation simulations invoke the `lib_size_effect_dynamics` method:
```python
exprU_O_L, exprS_O_L = sim.outlier_effect_dynamics(exprU_O, exprS_O, mean, scale)
```

3. **Dropouts**: 

In steady-state simulations invoke the `dropout_indicator` method:
```python
binary_ind = sim.dropout_indicator(expr_O_L, shape, percentile)
expr_O_L_D = np.multiply(binary_ind, expr_O_L)
```

In differentiation simulations invoke the `dropout_indicator_dynamics` method:
```python
binary_indU, binary_indS = sim.dropout_indicator_dynamics(exprU_O_L, exprS_O_L, shape, percentile)
exprU_O_L_D = np.multiply(binary_indU, exprU_O_L)
exprS_O_L_D = np.multiply(binary_indS, exprS_O_L)
```

4. **mRNA Count Matrix**: 

In steady-state simulations invoke the `convert_to_UMIcounts` method:
```python
count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)
```

In differentiation simulations invoke the `convert_to_UMIcounts_dynamics` method:
```python
count_matrix_U = sim.convert_to_UMIcounts_dynamics(exprU_O_L_D)
count_matrix_S = sim.convert_to_UMIcounts_dynamics(exprS_O_L_D)
```

The output of each of these modules including the "count matrix conversion" module are 3d numpy arrays of size (#cell_types * #gene * #cells_per_type). To convert them into a 2d expression matrix invoke numpy.concatenate as shown before. 

## Repository Contents
* SERGIO/ contains the python codes required for simulations.

* data_sets/ cotains 11 data sets including 6 steady-state and 5 differentiation simulated data. Each data set's folder contains the input files used in simulations as well the ground truth (gt) GRN. Differentiation data sets' folders also contain the differentiation graph (bMat) used in simulations.

* GNW_sampled_GRNs/ contains four networks sampled from the known regulatory network in Ecoli and Yeast using GeneNetWeaver (doi: 10.1093/bioinformatics/btr373). These networks might contain auto-regulatory edges and cycles.

* Demo/ contains demo input files for both steady-state and differentiation simulations. It also contains a jupyter notebook that runs demo simulations. Expected run time on a normal desktop computer for demo steady-state simulation is about 150 seconds and for demo differentiation simulations is about 120 seconds. 
