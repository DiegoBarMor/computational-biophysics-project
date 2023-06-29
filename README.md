# Computational Biophysics Project
This is a project of Computational Biophysics course 2022-2023 (University of Trento). The repository contains the source code and initial input files necessary to reproduce all the steps detailed by the [report](https://diegobarmor.github.io/computational-biophysics-project/).

The focus of this project is to analyse the dynamics of a protease of interest. The software used to simulate the trajectory is <code>GROMACS</code>, which integrates the equations of motion by means of the Verlet algorithm. Once the trajectories were obtained they were analysed, focusing on general aspects of the trajectories, as well as frames of interest, the active sites and the distribution of water around the protein. The analysis is mostly performed via Python 3 scripts. The library <code>MDAnalysis</code> facilitates dealing with GRO and XTC files, as well as performing different analysis methodologies.

![Protein trajectory](docs/animations/protein_opening.gif)
