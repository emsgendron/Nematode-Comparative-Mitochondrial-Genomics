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
	
	### Extract taxa from the database
		* Insert your database path
  		* Change taxa names to extract taxa of interest
  		* Insert your gene list path
	`while read -r line
	do
	geneseqs=$(grep --no-group-separator -A 1 -w "$line" $PATH_to_your_database) 
	grep --no-group-separator -A 1 $taxa <(echo "$geneseqs") >$taxa_"$line".fasta
	done<$PATH_to_genelist`
	
	
	### Filter fasta files to only include unique sequence IDs
		* change for directory and file name format for your fasta files for each gene
	`for file in *.fasta`
	`do`
	`grep '>' "$file" | sort | uniq >uniq_ids.list`
	`while read -r line`
	`do`
	`grep --no-group-separator -A 1 -w "$line" "$file" | head -n2`
	`done<uniq_ids.list >uniq_"$file"`
	`rm uniq_ids.list`
	`done`
	
	### Identify lowest number of shared taxa IDs across sequence files
		* We use COX1 as my starting point since it has the best coverage and then filter down from there as I work through each gene
		* Change inputs as you iterate through each gene
	`mkdir fastas_for_concatanation
	grep '>' uniq_*COX1.fasta | awk -F"|" '{print $1}' | sed 's/>//' | sort | uniq >fastas_for_concatanation/cox1_seqids.list 
	for file in uniq_*.fasta
	do
	ids=$(cat fastas_for_concatanation/cox1_seqids.list)
	while read -r line
	do
	grep -w "$line" $file | head -n1
	done<fastas_for_concatanation/cox1_seqids.list | sort | uniq | wc -l
	echo "$file"
	done
	grep --no-group-separator -f fastas_for_concatanation/cox1_seqids.list uniq_*$GENE.fasta  | awk -F"|" '{print $1}' | sed 's/>//' | sort | uniq >fastas_for_concatanation/lowest_taxaids.list3`
	
	
	### Subset sequence files once shared unique sequences have been identified
	`for file in *.fasta
	do
	ids=$(cat fastas_for_concatanation/lowest_taxaids.list2 | awk -F"|" '{print $1}' | sed 's/>//' | sort | uniq)
	while read -r line
	do
	grep -w -A 1 "$line" "$file" | head -n2
	done< <(echo "$ids") >fastas_for_concatanation/sub_"$file"
	done`
	
	### Check that all files have sequences in the same order and place across all gene files
	`for file in sub_*; do head -n73 "$file" | tail -n1; done`



	
