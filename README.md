# SERGIO v2 (Single-cell ExpRession of Genes In silicO)

Saurabh Sinha’s Lab, University of Illinois at Urbana-Champaign [Sinha Lab](https://www.sinhalab.net/sinha-s-home)

Developed by Payam Dibaeinia

## Description
SERGIO is a simulator for single-cell expression data guided by gene regulatory networks. SERGIO v2 is up to 100X faster than the v1 version. SERGIO v2 essentially simulates the same stochastic differential equations as v1 while provides users with additional functionalities.


## Getting Started
Install SERGIO from PyPI using:

```pip install sergio-scSIM```


## Usage

To run simulations a GRN and master regulators' profile are required.

### Step1: Build a GRN

A GRN is defined by a list of interactions. Each interactions is defined by a regulator and a target gene and by hill function's parameters k, n, and h. Interactions could be written in a .csv file with the following 9 columns:

* index: interactions index/id
* reg: regulator gene name
* coop: where the interactions is cooperative (1; includes more than one regulator) or non-cooperative (0; includes only one regulator)
* tar: target gene name
* k: interaction strength (use +/- values for activation/repression)
* n: non-linearity of the hill function
* h: half response
* reg_decay: decay rate of the regulator gene
* tar_decay: decay rate of the target gene

If you have a such a .csv file, build GRN using (make sure that your GRN does not have cycles):

```python
from SERGIO.GRN import grn_from_file
grn = grn_from_file(path = "your-file.csv")
```

Alternatively, if your GRN is not parametrized and you only have a .csv with the first four columns (index, reg, coop, tar) or if you want to re-parametrize your GRN, use:

```python
from SERGIO.GRN import grn_from_file
grn = grn_from_file(path = "your-file.csv", parametrize = True)
```

Alternatively, if you do not have any GRN, you can sample it (with certain number of genes) from a curated human network using:

```python
from SERGIO.GRN import grn_from_human
grn = grn_from_human(nGene = number-of-genes)
```

or you could sample it from Ecoli network:

```python
from SERGIO.GRN import grn_from_Ecoli
grn = grn_from_Ecoli(nGene = number-of-genes)
```

Finally, if you are willing to convert an interaction file from SERGIO v1 to a GRN for v2, use:

```python
from SERGIO.GRN import grn_from_v1
grn = grn_from_v1(path = "v1-interaction-file")
```

### Step2: Build MR profile

Once the GRN is built, its master regulators (MR) (i.e. the regulators that are not regulated by any up-stream regulator) can be extracted. Each cell type, is defined by the production rate of MRs in that cell type. Therefore, for simulating a certain number of cell types, the profile (production rates) of MRs in each cell type should be known. Use the following command to initialize a MR profile:

```python
from SERGIO.MR import mrProfile
mrs = grn.get_mrs()
mr_profs = mrProfile(MR_names = mrs, n_types = number-cell-types-to-simulate)
```

Next, you need to populate this profile with production rates of all MRs. If you have a .csv file containing the production rates of MRs (rows) in cell types (columns), you can use:

```python
mr_profs.build_from_file(path = "your-file.csv", header = 0)
```

Alternatively, you can define various levels for production rate of each MR in a cell type (e.g. 'L': lowly expressed, 'H': highly expressed), and define a reasonable range for each level (e.g. 'L':[1, 2.5], 'H': [3.5,5]). You can randomly assign a level to each MR in each cell types and randomly sample a production rate from the corresponding range. This process can be done with the following command

```python
mr_profs.build_rnd(range_dict={'L': [1, 2.5], 'H': [3.5, 5]})
```

Alternatively, if you wish to have certain number of markers per each cell type and/or certain number of MRs being highly expressed in all cell types, or certain number of MRs being lowly expressed in all cell types, use:

```python
mr_profs.build_rnd_complex(nMarker_per_type, nLow, nHigh)
```

### Step3: Simulate

Now, you are ready to start simulations (Note that if you used "grn_from_file" with, "parametrize = False", you have to set "update_half_resp = True" below):

```python
from SERGIO import sergio
grn.init(mr_profs, update_half_resp = True)
sim = sergio(grn)
sim.simulate(nCells = 200, noise_s = 1, safety_iter = 150, scale_iter = 10)
```
