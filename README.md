# Nematode-Comparative-Mitochondrial-Genomics
Pipeline for sorting, comparing, and extracting information on nematode mitochondrial genomes.

# REQUIRED SOFTWARE
* BASH compatable operating system or emulator
	* Linux
	* MACOSX
  	* Unix emulator (https://www.puttygen.com/windows-terminal-emulators)

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
		
		
		# This will require the installation of Miniconda/Anaconda for environment creation and activation
		https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

* PAML
http://abacus.gene.ucl.ac.uk/software/paml.html

* IQ-TREE # Can be either installed locally or you can use their webserver
http://www.iqtree.org/
# This pipeline makes use of the iq-tree webserver

* R Statistics and Plotting
https://www.r-project.org/
	# Required libraries
	'vegan'
	'ape'
	'tidyverse'
	'ggplot2'
	'viridis'
	'geiger'
	'pheatmap'
	'RColorBrewer'
	'lattice'

