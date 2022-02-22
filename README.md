## Maximum Entropy Model for Enzyme (MEME)

The code was initially developed by W.J Xie and B. Zhang at MIT studying the statistical physics of epigenome ([W.J. Xie and B. Zhang, Biophys Journal, 2019, 116, 2047-2056](https://www.sciencedirect.com/science/article/pii/S0006349519303066)). In the current version, the code was slightly modified to study protein evolution ([W.J. Xie, M. Asadi, and A. Warshel, PNAS, 2022, 119, e2122355119](https://www.pnas.org/content/119/7/e2122355119)).

### Descriptions
Programs to parameterize the MaxEnt model. Three techniques are used to accelerate the optimization: replica-exchange MCMC, MPI, and momentum-assisted SGD optimizer.

### Dependencies
GFortran v8.3.0, MPICH v3.3.2

### Input
statistics from MSA (experimental_constraints.txt; single-body frequencies are followed by double-body frequencies.)

### Output
parameters of the MaxEnt model (IsingHamiltonian_field.txt, IsingHamiltonian_coupling.txt in the params subfolder), the MaxEnt energy for mutated sequences (msa_mut_MaxEnt_energy.txt)
 
### Modules
init.f90 (initialize parameters, in the params subfolder), global.f90 (global variables used in all modules), hamiltonian.f90 (calculate configuration energies), mc_sampling.f90 (MC simulation module), ensemble_average.f90 (module to calculate statistics from simulated MSA), main.f90 (main module for the parameterization), energy.f90 (MaxEnt energy calculation given the mutated sequences). Due to size limit, the intermediate results for a specific example are not provided.

### How to use
'gfortran init.f90' (in the params subfolder. It can also further optimize parameters obtained from PLL approximation); 'bash sub.sh'; './energy'
