# Nematode-Comparative-Mitochondrial-Genomics
Pipeline for sorting, comparing, and extracting information on nematode mitochondrial genomes.

## REQUIRED SOFTWARE
* BASH compatable operating system or emulator
	* Linux
	* MACOSX
  	* Unix emulator (https://www.puttygen.com/windows-terminal-emulators)
* Miniconda/Anaconda
	* Software link: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
* PartitionFinder2
	* Software link: http://www.robertlanfear.com/partitionfinder/
	* This will require creating and activating a python2.7 environment for ParitionFinder
		* `conda create -n partitionfinder Python=2.7`
  		* `conda activate partitionfinder` 
	* Also install the following python packages:
		* conda, install, numpy, pandas, pyparsing, and scipy	
		* Example installation:
   			* `pip install -U scikit-learn`
			* `conda install pytables`

* PAML
	* Software link: http://abacus.gene.ucl.ac.uk/software/paml.html

* IQ-TREE
	* Software link: http://www.iqtree.org/
	* Can be either installed locally or you can use their webserver (this pipeline uses the webserver)
	

* R Statistics and Plotting
	* Software link: https://www.r-project.org/
	* Required libraries: vegan, ape, tidyverse, ggplot2, viridis, geiger, pheatmap, RColorBrewer, and lattice

 ## AGGREGATE AND CHECK GENOMES
* Obtain Gene Sequences from the mtMG-Database
* Download database fasta file: https://github.com/WormsEtAl/Nematode-Mitochondrial-Database
* Generate list of genes you are targeting and save as genelist
	* For example, the suggested list of genes for mitochondrial sequences: COX1, COX2, COX3, NAD1, NAD2, NAD3, NAD4, NAD4L, NAD5, NAD6, CYTB, ATP6
	
