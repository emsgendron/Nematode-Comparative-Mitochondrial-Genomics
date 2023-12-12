# Nematode-Comparative-Mitochondrial-Genomics
Pipeline for sorting, comparing, and extracting information on nematode mitochondrial genomes.

# REQUIRED SOFTWARE
* BASH compatable operating system or emulator
	* Linux
	* MACOSX
  	* Unix emulator (https://www.puttygen.com/windows-terminal-emulators)
* Miniconda/Anaconda
	* Software link: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
* PartitionFinder2
	* This will require creating and activating a python2.7 environment for ParitionFinder
		* `conda create -n partitionfinder Python=2.7`
  		* `conda activate partitionfinder` 
	* Software link: http://www.robertlanfear.com/partitionfinder/
	* Also install the following python packages:
		* conda, install, numpy, pandas, pyparsing, and scipy	
		* Example installation:
   			* `pip install -U scikit-learn`
			* `conda install pytables`

* PAML
	* Software link: http://abacus.gene.ucl.ac.uk/software/paml.html

* IQ-TREE
	* Can be either installed locally or you can use their webserver (this pipeline uses the webserver)
	* Software link: http://www.iqtree.org/

* R Statistics and Plotting
	* Software link: https://www.r-project.org/
	* Required libraries: vegan, ape, tidyverse, ggplot2, viridis, geiger, pheatmap, RColorBrewer, and lattice
	
