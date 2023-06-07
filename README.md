Snakemake Workflow for Collisional Cross Section Predictions from SMILES Input

This repository contains a Snakemake workflow manager for predicting collisional cross sections (CCS) using simplified molecular-input line-entry system (SMILES) structures. The workflow allows users to predict CCS values for multiple protonated/deprotonated adducts and models with high automation and parallelized computation on high-performance computing (HPC) systems.

Citation: If you use this workflow in your research, please consider citing the following publication:

Susanta Das, Kiyoto Aramis Tanemura, Laleh Dinpazhoh, Mithony Keng, Christina Schumm, Lydia Leahy, Carter K Asef, Markace Rainey, Arthur S. Edison, Facundo M. Fernández, and Kenneth M. Merz Jr. "In Silico Collision Cross Section Calculations to Aid Metabolite Annotation". J. Am. Soc. Mass Spectrom. 2022, 33, 5, 750–759. [DOI: 10.1021/jasms.1c00315](https://doi.org/10.1021/jasms.1c00315)

Installation Instructions:
Before running the Snakemake workflow, ensure that you have the following prerequisite software installed:
Dimorphite-DL: For ionization state determination. Available at: [https://durrantlab.pitt.edu/dimorphite-dl](https://durrantlab.pitt.edu/dimorphite-dl)
RDKit: For conformation generation. Available at: [https://www.rdkit.org](https://www.rdkit.org)
ASE-ANI: For conformation filtering. Available at: [https://github.com/isayev/ASE_ANI](https://github.com/isayev/ASE_ANI)
QUICK: For quantum calculations. Available at: [https://github.com/merzlab/QUICK](https://github.com/merzlab/QUICK)
hpccs: For CCS calculation. Available at: [https://github.com/cepid-cces/hpccs](https://github.com/cepid-cces/hpccs)

To install Snakemake, follow these steps:
1. Activate your base conda environment:
```
conda activate base
```
2. Create a new conda environment named "snakemake" and install Snakemake:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
3. Activate the "snakemake" environment:
```
conda activate snakemake
```
Usage Guide:
1. Specify the paths to the installed software in the `paths.json` file.
2. Configure the walltime, memory requirement, and node specifications in the `cluster.yaml` file.
3. Provide all the necessary arguments in the `arguments.json` file.
4. To run the workflow, execute the following command:
```
bash scheduler.sh
```
5. After successful job completion, the workflow will generate a `ccs.txt` file in the `ensemble` (or `ensemble_fast`) folder. The `ccs.txt` file contains the Boltzmann average CCS values computed using this workflow.

Web server based CCS Calculation:
In addition to running the Snakemake workflow locally, you can also access a web-based CCS calculation tool provided by [POMICS](https://www.pomics.org/). This website offers free academic usage and calculates CCS upon user request. The website runs the snakemake_ccs workflow at the backend and provides CCS results along with all the metadata.
Please visit POMICS to access the CCS calculation tool and explore further capabilities.

License:
This Snakemake workflow for collisional cross section predictions is licensed under an academic license granted by Michigan State University. By downloading the software through the GitHub repository, you agree to the terms of the academic license. For more information about the snakemake_ccs tool, please contact Kenneth M. Merz, Jr. at merz@msu.edu or Susanta Das at dassusan@msu.edu.


