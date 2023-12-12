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
		while read -r line
		do
		geneseqs=$(grep --no-group-separator -A 1 -w "$line" $PATH_to_your_database)  
		grep --no-group-separator -A 1 $taxa <(echo "$geneseqs") >$taxa_"$line".fasta 
		done<$PATH_to_genelist
	* Insert your database path
  	* Change taxa names to extract taxa of interest
  	* Insert your gene list path
	
	
	### Filter fasta files to only include unique sequence IDs
		for file in *.fasta
		do
		grep '>' "$file" | sort | uniq >uniq_ids.list
			while read -r line
			do
			grep --no-group-separator -A 1 -w "$line" "$file" | head -n2
			done<uniq_ids.list >uniq_"$file"
		rm uniq_ids.list
		done
  	* Change directory and file extension format to match the location and names for your sequence files extracted above.
	
	### Identify lowest number of shared taxa IDs across sequence files
	* We use COX1 as my starting point since it has the best coverage and then filter down from there as I work through each gene
  	###
   		mkdir fastas_for_concatanation

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
	
 	* Change inputs as you iterate through each gene
	###
  		grep --no-group-separator -f fastas_for_concatanation/cox1_seqids.list uniq_*$GENE.fasta  | awk -F"|" '{print $1}' | sed 's/>//' | sort | uniq >fastas_for_concatanation/lowest_taxaids.list3`
	
	
	### Subset sequence files once shared unique sequences have been identified
		for file in *.fasta
		do
		ids=$(cat fastas_for_concatanation/lowest_taxaids.list2 | awk -F"|" '{print $1}' | sed 's/>//' | sort | uniq)
			while read -r line
			do
			grep -w -A 1 "$line" "$file" | head -n2
			done< <(echo "$ids") >fastas_for_concatanation/sub_"$file"
		done
	
	### Check that all files have sequences in the same order and place across all gene files
	`for file in sub_*; do head -n73 "$file" | tail -n1; done`


## FORMAT SEQUENCE IDS FOR SHORTHAND
	mkdir downstream_formatted
	cp sub* downstream_formatted
	cd downstream_formatted

	for file in sub*
	do
	seqids=$(grep '>' "$file")
	numseqs=$(wc -l <(echo "$seqids") | cut -f1 -d" ")
		fmtseqids=$(for i in $( eval echo {1..$numseqs})
		do
		echo Seq_"$i"
		done)
	echo "$seqids" >seqids.tmp
	echo "$fmtseqids" >fmtseqids.tmp
	paste seqids.tmp fmtseqids.tmp >"$file"_ID_map.tsv
		while read -r line;
		do
		toreplace=$(awk -F'\t' '{print $1}' <(echo "$line"))
		replacewith=$(awk -F'\t' '{print $2}' <(echo "$line"))
		orgID=$(grep "$toreplace" "$file")
		echo "$orgID"
		orgIDmap=$(grep -w "$replacewith" "$file"_ID_map.tsv | awk -F'\t' '{print $1}')
		echo "$orgIDmap"
		newIDmap=$(grep -w "$replacewith" "$file"_ID_map.tsv | awk -F'\t' '{print $2}')
		echo "$newIDmap"
		newID=$(sed "s/$toreplace/$replacewith/" "$file" | grep -w "$replacewith")
		echo "$newID"
		if [ "$orgID" != "$orgIDmap" ]; then
		echo "Original Map IDs do not match file's original sequence IDs"
		echo $file
		exit 0
		else
		echo "Original IDs Match"
		fi
		if [ "$newID" != "$newIDmap" ]; then
		echo "New Map IDs do not match file's updated sequence IDs"
		echo $file
		exit 0
		else
		echo "New IDs Match"
		fi
		sed -i "s/$toreplace/>$replacewith/" "$file"
		done<"$file"_ID_map.tsv
	done
	
 	rm seqids.tmp
	rm fmtseqids.tmp

## ALIGN SEQUENCES
	for file in sub*.fasta
	do
	mafft --auto --adjustdirectionaccurately "$file" >adjacc_mafft_"$file"
	done

## CONCATONATE ALIGNMENTS TO GENERATE "GENOME" SEQUENCE
### Nucleotides "genome"
* Ensure sequences are formated as two-line interleaved
* format sequence IDs for any potentially rearranged sequences
##
	for file in adjacc_mafft_*
	do
	sed -i 's/_R_//g' "$file"
	done

* Concatenate Sequence Files
##
	paste adjacc*.fasta | sed 's/\t>Seq.*//g' | sed 's/\t//g' >concatanated_aln_pcgs.fasta
	grep 'Seq' concatanated_aln_pcgs.fasta | sort | uniq >uniq_concat_ids.list

* Insure no duplicate sequence IDs are present
##
	while read -r line
	do
	grep --no-group-separator -A 1 -w "$line" concatanated_aln_pcgs.fasta | head -n2
	done<uniq_concat_ids.list >uniq_concatanated_aln_pcgs.fasta

### Amino Acid "genome"
* Insure sequences are formated as two-line interleaved
* Concatenate Sequence Files
##
	paste *.pep.aln | sed 's/\t>Seq.*//g' | sed 's/\t//g' >concatanated_aln_pcgs.faa
	grep 'Seq' concatanated_aln_pcgs.faa | sort | uniq >uniq_concat_ids.list

* Insure no duplicate sequence IDs are present
##
	while read -r line
	do
	grep --no-group-separator -A 1 -w "$line" concatanated_aln_pcgs.faa | head -n2
	done<uniq_concat_ids.list >uniq_concatanated_aln_pcgs.faa

* Move your concatenated genomes into their own directory
##
	mkdir concat_genomes
	mv concatanated* concat_genomes/
	mv uniq_concat* concat_genomes/



	
