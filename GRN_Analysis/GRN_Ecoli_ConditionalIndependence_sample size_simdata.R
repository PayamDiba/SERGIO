############################################################################
# Conditional independencies of a GRN vs sample size using simulated data  # 
############################################################################

# Load packages
library(dagitty)
library(lavaan)
library(sna)
library(DOT)
library(igraph)
library(pcalg)
library(SEMgraph)
library(reticulate)
library(dplyr)


#########################
# Sergio: python from R #
#########################

setwd("C:/Users/moha456/OneDrive - PNNL/Documents/Missing data//SERGIO/Demo")
use_python("C:/Users/moha456/Anaconda3")
os <- import("os")
os$listdir(".")

# Sourcing python scripts
source_python("gene.py")
source_python("sergio.py")

# Import modules
np <- import("numpy", convert = FALSE)
pd <- import("pandas", convert = FALSE)
sg <- import("sergio", convert = FALSE)

# Simulate Clean Data _ Steady State Simulation for single cell data
sim = sg$sergio(number_genes= as.integer(100), number_bins = as.integer(9), number_sc = as.integer(300),
                noise_params = 1, sampling_state=15, decays= 0.8, noise_type='dpd')
sim$build_graph(input_file_taregts ='steady-state_input_GRN.txt', input_file_regs='steady-state_input_MRs.txt', shared_coop_state=2)
sim$simulate()
expr = sim$getExpressions()
expr_clean_ss = np$concatenate(expr, axis = as.integer(1))
head(expr_clean_ss)

# Simulate Clean Data _ Steady State Simulation for single cell data
sim = sg$sergio(number_genes= as.integer(100), number_bins = as.integer(1), number_sc = as.integer(300),
                noise_params = 1, sampling_state=15, decays= 0.8, noise_type='dpd')
sim$build_graph(input_file_taregts ='single_cell_GRN.txt', input_file_regs='single_cell_MRs.txt', shared_coop_state=2)
sim$simulate()
expr = sim$getExpressions()
expr_clean_ss = np$concatenate(expr, axis = as.integer(1))
head(expr_clean_ss)
