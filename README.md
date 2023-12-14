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

* MUSCLE Sequence Aligner
	* Software link: https://www.drive5.com/muscle/  

* R Statistics and Plotting
	* Software link: https://www.r-project.org/
	* Required libraries: vegan, ape, tidyverse, ggplot2, viridis, geiger, pheatmap, RColorBrewer, FSA, and lattice

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
	for file in adjacc_mafft_*
	do
	sed -i 's/_R_//g' "$file"
	done
* Ensure sequences are formated as two-line interleaved
* format sequence IDs for any potentially rearranged sequences
#
	paste adjacc*.fasta | sed 's/\t>Seq.*//g' | sed 's/\t//g' >concatanated_aln_pcgs.fasta
	grep 'Seq' concatanated_aln_pcgs.fasta | sort | uniq >uniq_concat_ids.list
* Concatenate Sequence Files
#
	while read -r line
	do
	grep --no-group-separator -A 1 -w "$line" concatanated_aln_pcgs.fasta | head -n2
	done<uniq_concat_ids.list >uniq_concatanated_aln_pcgs.fasta
* Insure no duplicate sequence IDs are present
### Amino Acid "genome"
	paste *.pep.aln | sed 's/\t>Seq.*//g' | sed 's/\t//g' >concatanated_aln_pcgs.faa
	grep 'Seq' concatanated_aln_pcgs.faa | sort | uniq >uniq_concat_ids.list
* Insure sequences are formated as two-line interleaved
* Concatenate Sequence Files
#
	while read -r line
	do
	grep --no-group-separator -A 1 -w "$line" concatanated_aln_pcgs.faa | head -n2
	done<uniq_concat_ids.list >uniq_concatanated_aln_pcgs.faa
* Insure no duplicate sequence IDs are present
#
	mkdir concat_genomes
	mv concatanated* concat_genomes/
	mv uniq_concat* concat_genomes/
* Move your concatenated genomes into their own directory


## CODON ALIGN EACH GENE FASTA
* Now we do the codon alignment after determining correct orientation for each sequence
* remove alignment gaps and realign as a codon alignment

### Format mafft alignment files for CODON alignment
	for file in mafftadj_*
	do
	sed -i 's/-//g' "$file"
	done
* SKIP IF NOT APPLYING TO CODON ALIGNMENT

### Codon alignment
	for file in mafftadj_*
	do
	perl ~/perl_scripts/codon_alignment-master/codon_alignment_e1.pl "$file" 5 0.1
	wait
	muscle -align "$file".pep -output "$file".pep.aln
	wait
	perl ~/perl_scripts/codon_alignment-master/codon_alignment_e1.pl "$file" 5 0.1
	wait
	muscle -align "$file".pep.aln.cleanup -output "$file".pep.aln.cleanup.aln
	wait
	perl ~/perl_scripts/codon_alignment-master/codon_alignment_e1.pl "$file" 5 0.1
	wait
	done




## CONVERT SEQUENCE FILES TO PHYLIP FORMAT

### Create Nucleotide files
	for file in *.phylip
	do
	seq=$(head -n2 "$file" | tail -n1)
	awk '{printf "%s%s", $0, (/Seq/ ? " " : ORS)}' "$file" >"$file".phylip2
	sed -i 's/ /\t/' "$file".phylip2
	done
	rm *.tmp

### Create PAML formatted phylip files
	for file in *.phylip2
	do
	sed 's/ /  /' "$file" | sed 's/\t/  /' >"$file".paml
	done

### Create Amino acid files
	for file in *.aln
	do
	seqids=$(grep 'Seq' "$file")
	numseqs=$(wc -l <(echo "$seqids") | cut -f1 -d" ")
	seq=$(head -n2 "$file" | tail -n1)
	echo "$numseqs@${#seq}" >header.tmp
	cat header.tmp "$file" >"$file".phy.tmp
	sed -i 's/>//' "$file".phy.tmp
	sed -i 's/@/\t/' "$file".phy.tmp
	awk '{printf "%s%s", $0, (/Seq/ ? " " : ORS)}' "$file".phy.tmp >"$file".phy2
	sed -i 's/ /\t/' "$file".phy2
	done
	rm *.tmp



# CREATE PARTITION CONFIGURATION FILES 
## for nucleotide files
	for file in *.phylip
	do
	seqlength=$(head -n1 "$file" | awk -F'\t' '{print $2}')
	gene=$(echo $file | awk -F"_" '{print $NF}' | awk -F"." '{print $1}')
	output=$(echo alignment = "$file"";"
	echo branchlengths = linked";"
	echo models = all";"
	echo model_selection = AICc";"
	echo [data_blocks]
	echo "$gene"_pos1=1-"$seqlength""\3;"
	echo "$gene"_pos2=2-"$seqlength""\3;"
	echo "$gene"_pos3=3-"$seqlength""\3;"
	echo [schemes]
	echo search = greedy";")
	echo "$output" >"$file".cfg
	done

## for amino acid files
	for file in *.phy
	do
	seqlength=$(head -n1 "$file" | awk -F'\t' '{print $2}')
	gene=$(echo $file | awk -F"_" '{print $NF}' | awk -F"." '{print $1}')
	output=$(echo alignment = "$file"";"
	echo branchlengths = linked";"
	echo models = all";"
	echo model_selection = AICc";"
	echo [data_blocks]
	echo "$gene"_pos1=1-"$seqlength""\3;"
	echo "$gene"_pos2=2-"$seqlength""\3;"
	echo "$gene"_pos3=3-"$seqlength""\3;"
	echo [schemes]
	echo search = greedy";")
	echo "$output" >"$file".cfg2
	done

## Collect info for each gene for concatenated partition
	for file in *.phylip; do seq=$(tail -n1 "$file"); echo "${#seq}"; done

### EXAMPLE FORMAT
	alignment = concatanated_aln_pcgs.phylip.phylip2;
	branchlengths = linked; 
	models = all; 
	model_selection = AICc; 
	
	[data_blocks]
	ATP6_pos1=1-744\3;
	ATP6_pos2=2-744\3;
	ATP6_pos3=3-744\3;
	COX1_pos1=745-3789\3;
	COX1_pos2=746-3789\3;
	COX1_pos3=747-3789\3;
	COX2_pos1=3790-4935\3;
	COX2_pos2=3791-4935\3;
	COX2_pos3=3792-4935\3;
	COX3_pos1=4936-6402\3;
	COX3_pos2=4937-6402\3;
	COX3_pos3=4938-6402\3;
	CYTB_pos1=6403-8883\3;
	CYTB_pos2=6404-8883\3;
	CYTB_pos3=6405-8883\3;
	NAD1_pos1=8884-10479\3;
	NAD1_pos2=8885-10479\3;
	NAD1_pos3=8886-10479\3;
	NAD2_pos1=10480-12039\3;
	NAD2_pos2=10481-12039\3;
	NAD2_pos3=10482-12039\3;
	NAD3_pos1=12040-12624\3;
	NAD3_pos2=12041-12624\3;
	NAD3_pos3=12042-12624\3;
	NAD4_pos1=12625-14739\3;
	NAD4_pos2=12626-14739\3;
	NAD4_pos3=12627-14739\3;
	NAD4L_pos1=14740-15204\3;
	NAD4L_pos2=14741-15204\3;
	NAD4L_pos3=14742-15204\3;
	NAD5_pos1=15205-17529\3;
	NAD5_pos2=15206-17529\3;
	NAD5_pos3=15207-17529\3;
	NAD6_pos1=17530-18300\3;
	NAD6_pos2=17531-18300\3;
	NAD6_pos3=17532-18300\3;
	
	[schemes]
	search = greedy;

# IDENTIFY SEQUENCE PARTITIONS FOR RATE MODELS
## Activate a python2.7 specific environment
	conda activate partitionfinder
## format input phylip file names in the config file
	for file in *.phylip.cfg
	do
	sed -i .bak 's/phylip/phylip.phylip2/' "$file"
	rm *.bak
	done


### for Nucleotides
	mkdir fnas
	for file in *cfg
	do
	gene=$(echo $file | awk -F"_" '{print $NF}' | awk -F"." '{print $1}')
	mv "$file" partition_finder.cfg
	wait
	python2.7 $PATH_to_PartitionFinder/PartitionFinder.py $PATH_for_output
	mv analysis "$gene"_analysis
	mv "$gene"_analysis fnas/
	mv partition_finder.cfg "$file"
	done
* change directory as needed

### for Amino Acids
	mkdir faas
	for file in *.phy.cfg2
	do
	gene=$(echo $file | awk -F"_" '{print $NF}' | awk -F"." '{print $1}')
	mv "$file" partition_finder.cfg
	python2.7 $PATH_to_PartitionFinder/PartitionFinder.py $PATH_for_output
	mv analysis "$gene"_analysis
	mv "$gene"_analysis faas/
	mv partition_finder.cfg "$file"
	done
* change directory as needed

## Extract nexus formatted partition 
	for dir in *_analysis
	do
	gene=$(echo $dir | cut -f1 -d"_")
	grep -A 6 "#nexus" "$gene"_analysis/best_scheme.txt >"$gene"_partition.txt
	done


# PARSING PAML OUTPUTS: MUTATION RATES RELATIVE TO OUTGROUP
* change header length for the number of comparisons you want to extract
* CHANGE PARTITION OUTGROUP COMPARISON IDS BASED ON DESIRED OUTGROUP (e.g. Seq_99 is the ID of the sequence used as an outgroup for comparison
	- Each line sets the outgroups you wish to compare to and extracts line of mutation metrics where each outgroup sequence is present
##
	while read -r line2
	do
	header=$(echo "SeqID"@"Out01"@"Out02"@"Out3"@"Out4"@"Out5"@"Avg")
	sed 's/@/\t/g' <(echo "$header")
	while read -r line;
	do
	
	rates08=$(grep -A 6 \(Seq_99\) "$line2"/"$line2"_branch_w1_mlc | grep -A 6 -w \("$line" | tail -n1 | awk -F"dN/dS" '{print $2}')
	rates09=$(grep -A 6 \(Seq_100\) "$line2"/"$line2"_branch_w1_mlc | grep -A 6 -w \("$line" | tail -n1 | awk -F"dN/dS" '{print $2}')
	rates11=$(grep -A 6 \(Seq_101\) "$line2"/"$line2"_branch_w1_mlc | grep -A 6 -w \("$line" | tail -n1 | awk -F"dN/dS" '{print $2}')
	rates20=$(grep -A 6 \(Seq_102\) "$line2"/"$line2"_branch_w1_mlc | grep -A 6 -w \("$line" | tail -n1 | awk -F"dN/dS" '{print $2}')
	rates33=$(grep -A 6 \(Seq_103\) "$line2"/"$line2"_branch_w1_mlc | grep -A 6 -w \("$line" | tail -n1 | awk -F"dN/dS" '{print $2}')
	
	# DNDS ratios
	dnds08=$(awk -F" " '{print $2}' <(echo "$rates08"))
	dnds09=$(awk -F" " '{print $2}' <(echo "$rates09"))
	dnds11=$(awk -F" " '{print $2}' <(echo "$rates11"))
	dnds20=$(awk -F" " '{print $2}' <(echo "$rates20"))
	dnds33=$(awk -F" " '{print $2}' <(echo "$rates33"))
	
	# DN rates
	dn08=$(awk -F"dN =" '{print $2}' <(echo "$rates08") | awk -F" " '{print $1}')
	dn09=$(awk -F"dN =" '{print $2}' <(echo "$rates09") | awk -F" " '{print $1}')
	dn11=$(awk -F"dN =" '{print $2}' <(echo "$rates11") | awk -F" " '{print $1}')
	dn20=$(awk -F"dN =" '{print $2}' <(echo "$rates20") | awk -F" " '{print $1}')
	dn33=$(awk -F"dN =" '{print $2}' <(echo "$rates33") | awk -F" " '{print $1}')
	
	# DS rates
	ds08=$(awk -F"dS" '{print $2}' <(echo "$rates08") | sed 's/=//')
	ds09=$(awk -F"dS" '{print $2}' <(echo "$rates09") | sed 's/=//')
	ds11=$(awk -F"dS" '{print $2}' <(echo "$rates11") | sed 's/=//')
	ds20=$(awk -F"dS" '{print $2}' <(echo "$rates20") | sed 's/=//')
	ds33=$(awk -F"dS" '{print $2}' <(echo "$rates33") | sed 's/=//')
	
	#### CHANGE OUTPUT STRINGS FOR EACH DESIRED RATE ####
	out=$(echo "$line"@"$dn08"@"$dn09"@"$dn11"@"$dn20"@"$dn33")
	sed 's/@/\t/g' <(echo "$out")
	#### CHANGE FILE NAME FOR EACH DESIRED RATE OUTPUT ####
	done<$PATH_to_outgroup_IDs_list >"$line2"_DN_rates
	done<genelist # use the same style gene list as before



# RENAME RATE TABLES WITH TAXONOMIES
* REQUIRES MAP TABLE WITH SIMPLIFIED SEQ IDS ROW MATCHED TO DESIRED TAXA IDS ###
## CHANGE FILE NAMES TO MATCH TYPE OF RATE FILE YOU ARE ALTERING ####
	while read -r line2;
	do
	while read -r line
	do replacewith=$(grep -w "$line" rhabditina_ID_map.tsv | cut -f1) # change ID map file based on needs
	sed -i ".bak" "s/$line\t/$replacewith\t/g" "$line2"_DN_rates
	awk -F"|" '{print $1}' "$line2"_DN_rates >specieslist
	done<seqlist
	paste specieslist "$line2"_DN_rates | sed 's/\t/,/g' >"$line2"_Rtable.csv
	done<genelist

# PARSING IQTREE RATE TABLES
## Renaming short IDs with full taxonomy
	while read -r line
	do
	tax=$(awk -F"\t" '{print $1}' <(echo "$line"))
	ID=$(awk -F"\t" '{print $2}' <(echo "$line"))
	sed -i .bak "s/$ID:/$tax:/" subset_alignment.fasta.contree
	done<subset_ID_map.tsv
* Change input tree as needed 
* Change input seqID map as needed


## Extract and output rate matrix
	for dir in *ina/iqtree_outputs
	do
	ratemat=$(grep -A 5 'Rate matrix' "$dir"/*.iqtree | tail -n4 | sed 's/ * /@/g')
	grep Command "$dir"/*log
		while read -r line
		do
		awk -F'@' '{print $3}' <(echo "$line")
		awk -F'@' '{print $4}' <(echo "$line")
		awk -F'@' '{print $5}' <(echo "$line")
		awk -F'@' '{print $6}' <(echo "$line")
		done< <(echo "$ratemat")
	done


## Extract and output State Frequencies
	for dir in *ina/iqtree_outputs
	do
	grep Command "$dir"/*log
	statfq=$(grep -A 5 'State frequencies' "$dir"/*.iqtree | tail -n4 | sed 's/ * /@/g')
		while read -r line
		do
		awk -F'@' '{print $4}' <(echo "$line")
		done< <(echo "$statfq")
	done

## Extract and output Rate parameters
	for dir in *ina/iqtree_outputs
	do
	grep Command "$dir"/*log
	R=$(grep -A 7 'Rate parameter' "$dir"/*.iqtree | tail -n6 | sed 's/ * /@/g')
		while read -r line
		do
		awk -F'@' '{print $3}' <(echo "$line")
		done< <(echo "$R")
	done

## Extract and output Invariable sites and gamma shapes for alpha parameter
	for dir in *ina/iqtree_outputs
	do
	grep Command "$dir"/*log
	grep 'Proportion of invariable sites' "$dir"/*.iqtree | awk -F': ' '{print $2}'
	grep 'Gamma shape alpha' "$dir"/*.iqtree | awk -F': ' '{print $2}'
	cat0=$(grep -A 9 'Model of rate heterogeneity' "$dir"/*.iqtree | tail -n5 | sed 's/ * /@/g')
		while read -r line
		do
		awk -F'@' '{print $3}' <(echo "$line")
		awk -F'@' '{print $4}' <(echo "$line")
		done< <(echo "$cat0")
	done

# PLOT MUTATION RATE DATA
## Load R in your terminal
## Load Libraries
	library(vegan)
	library(ape)
	library(tidyverse)
	library(ggplot2)
	library(viridis)
	library(geiger)
	library(pheatmap)
	library(RColorBrewer)
	library(lattice)

## Master metadata table
	mandat<-read.csv("../phylosorted_spirurina_metadata.csv", header=T, row.names=1)

## List of prefixes matching file names in working directory
	genes<-c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6")

* Rates Rtables Requires the following fields 
	* One column for comparisons to each out group member
	* One column of full sequence IDs with taxonomies
	* One column of the most basic but unique tax ID that overlaps with rownames in the mandat object

### Gene Relative Scale Gradients Script
	for( i in 1:length(genes)) {
	name<-genes[i]
	final.tab<-matrix(,ncol=238)
	filename<-paste(name,"_Rtable.csv", sep="")
	dnds_data <- read.csv(filename, header=F)
		for( a in 1:nrow(mandat)){
		species<-rownames(mandat)[a]
			for( x in 1:nrow(dnds_data)) {
			id<-dnds_data$V1[x]
			id.pos<-grep(id, species)
			rates<-dnds_data[,3:7]
			sub.rates<-apply(rates,1,as.numeric)
			averages<-apply(sub.rates,2,mean)
			dnds.tmp<-cbind(dnds_data,averages)
				if ( identical(id.pos,integer(0)) == FALSE ) {
				c1<-mandat[rownames(mandat) == species,]
				c2<-dnds.tmp[dnds.tmp$V1 == id,]
				point<-cbind(mandat[rownames(mandat) == species,],dnds.tmp[dnds.tmp$V1 == id,])
				colnames(final.tab)<-colnames(point)
				final.tab<-rbind(final.tab,point)
				}
				write.table(final.tab, paste("phylosorted_",name,".csv", sep=""), sep=",")
			}
		}
		print(paste("Gene_",i,":",name,"-Complete", sep=""))
	}
* WARNING: ncols(empty matrix) = ncol(master table) + ncol(rate table) + 1 (for adding the average rate values) + 1

### Fixed Scale Gradients Script
	for( i in 1:length(genes)) {
	name<-genes[i]
	final.tab<-matrix(,ncol=265)
	filename<-paste(name,"_Rtable.csv", sep="")
	dnds_data <- read.csv(filename, header=F)
		for( a in 1:nrow(mandat)){
		species<-rownames(mandat)[a]
			for( x in 1:nrow(dnds_data)) {
			id<-dnds_data$V1[x]
			id.pos<-grep(id, species)
			rates<-dnds_data[,3:6]
			sub.rates<-apply(rates,1,as.numeric)
			averages<-apply(sub.rates,2,mean)
			dnds.tmp<-cbind(dnds_data,averages)
				if ( identical(id.pos,integer(0)) == FALSE ) {
				c1<-mandat[rownames(mandat) == species,]
				c2<-dnds.tmp[dnds.tmp$V1 == id,]
				point<-cbind(mandat[rownames(mandat) == species,],dnds.tmp[dnds.tmp$V1 == id,])
				colnames(final.tab)<-colnames(point)
				final.tab<-rbind(final.tab,point)
				}
			}
			final.tab$V1 <- factor(final.tab$V1, levels = final.tab$V1)
			pdfname<-paste(name,"_averages.pdf",sep="")
			pdf(paste(pdfname))
			print(ggplot(final.tab, aes(x = 1, y = factor(V1), fill = as.numeric(averages))) + geom_tile() + scale_fill_viridis(limits=c(0,150)) + theme(text = element_text(size = 10)) + scale_y_discrete(limits=rev, position = 'right') + labs(fill="Average dN Relative to Outgroup"))
			dev.off()
		}
		print(max(final.tab$averages[2:length(final.tab$averages)]))
		print(paste("Gene_",i,":",name,"-Complete", sep=""))
	}
* Alter limits scale based on maximums in dataset
* WARNING: ncols(empty matrix) = ncol(master table) + ncol(rate table) + 1 (for adding the average rate values) + 1

## Plot from Full Metadata Table
	filename<-"phylosorted_table.csv"
	dat<-read.csv("phylosorted_table.csv", header=T)
	
	max(dat$DN.Rates, na.rm = TRUE)
	max(dat$DS.Rates, na.rm = TRUE)
	max(dat$DNDS.Rates, na.rm = TRUE)
	min(dat$DN.Rates, na.rm = TRUE)
	min(dat$DS.Rates, na.rm = TRUE)
	min(dat$DNDS.Rates, na.rm = TRUE)
	
	enoplea.dat<-subset(dat, dat$Class == "Enoplea")
	plot.dat<-apply(enoplea.dat[,3:5],c(1,2),as.numeric)
	plot.dat<-as.data.frame(plot.dat)
	plot.dat$Gene<-enoplea.dat$Gene
	plot.dat$Species<-enoplea.dat$Species
	factor.set<-plot.dat[plot.dat$Gene == "COX1",]
	plot.dat$Species <- factor(plot.dat$Species, levels = factor.set$Species)
	
	pdf("dn_enoplea.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DN.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,20)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN Relative to Outgroup")
	dev.off()
	
	pdf("ds_enoplea.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,310)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dS Relative to Outgroup")
	dev.off()
	
	pdf("dnds_enoplea.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DNDS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,3)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN:dS Relative to Outgroup")
	dev.off()
	
	tylenchina.dat<-subset(dat, dat$Suborder == "Tylenchina")
	plot.dat<-apply(tylenchina.dat[,3:5],c(1,2),as.numeric)
	plot.dat<-as.data.frame(plot.dat)
	plot.dat$Gene<-tylenchina.dat$Gene
	plot.dat$Species<-tylenchina.dat$Species
	factor.set<-plot.dat[plot.dat$Gene == "COX1",]
	plot.dat$Species <- factor(plot.dat$Species, levels = factor.set$Species)
	
	pdf("dn_tylenchina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DN.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,20)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN Relative to Outgroup")
	dev.off()
	
	pdf("ds_tylenchina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,310)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dS Relative to Outgroup")
	dev.off()
	
	pdf("dnds_tylenchina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DNDS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,3)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN:dS Relative to Outgroup")
	dev.off()
	
	spirurina.dat<-subset(dat, dat$Suborder == "Spirurina")
	plot.dat<-apply(spirurina.dat[,3:5],c(1,2),as.numeric)
	plot.dat<-as.data.frame(plot.dat)
	plot.dat$Gene<-spirurina.dat$Gene
	plot.dat$Species<-spirurina.dat$Species
	factor.set<-plot.dat[plot.dat$Gene == "COX1",]
	plot.dat$Species <- factor(plot.dat$Species, levels = factor.set$Species)
	
	pdf("dn_spirurina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DN.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,20)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN Relative to Outgroup")
	dev.off()
	
	pdf("ds_spirurina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,310)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dS Relative to Outgroup")
	dev.off()
	
	pdf("dnds_spirurina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DNDS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,3)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN:dS Relative to Outgroup")
	dev.off()
	
	rhabditina.dat<-subset(dat, dat$Suborder == "Rhabditina")
	plot.dat<-apply(rhabditina.dat[,3:5],c(1,2),as.numeric)
	plot.dat<-as.data.frame(plot.dat)
	plot.dat$Gene<-rhabditina.dat$Gene
	plot.dat$Species<-rhabditina.dat$Species
	factor.set<-plot.dat[plot.dat$Gene == "COX1",]
	plot.dat$Species <- factor(plot.dat$Species, levels = factor.set$Species)
	
	pdf("dn_rhabditina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DN.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,20)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN Relative to Outgroup")
	dev.off()
	
	pdf("ds_rhabditina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,310)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dS Relative to Outgroup")
	dev.off()
	
	pdf("dnds_rhabditina.pdf")
	ggplot(plot.dat, aes(x = Gene, y = as.factor(Species), fill = DNDS.Rates)) + 
	geom_tile() + 
	scale_fill_viridis(limits=c(0,3)) + 
	theme(text = element_text(size = 8), axis.text.x = element_text(angle = 90)) +
	scale_y_discrete(limits=rev, position = 'right') +
	labs(fill="Average dN:dS Relative to Outgroup")
	dev.off()

* Required reformatted table with single column variables
* Hand sorted taxa based on predicted phylogeny



# STATS
## Load required libraries
	library(FSA)

## import data
	dat<-read.csv("tablev3.csv", header=T)
* a table containing genome characteristics, skews, and substitution rates paired to each genome

## Taxonomy explaining genome characteristics
### GC Skew
	kruskal.test(dat$Genome.GC.Skew ~ Class, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Subclass, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Order, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Suborder, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Superfamily, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Family, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Genus, data = dat)
### Size
	kruskal.test(dat$Length..bp. ~ Class, data = dat)
	kruskal.test(dat$Length..bp. ~ Subclass, data = dat)
	kruskal.test(dat$Length..bp. ~ Order, data = dat)
	kruskal.test(dat$Length..bp. ~ Suborder, data = dat)
	kruskal.test(dat$Length..bp. ~ Superfamily, data = dat)
	kruskal.test(dat$Length..bp. ~ Family, data = dat)
	kruskal.test(dat$Length..bp. ~ Genus, data = dat)
### PCG GC Skew
	kruskal.test(dat$PCG.GC.Skew ~ Class, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Subclass, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Order, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Suborder, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Superfamily, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Family, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Genus, data = dat)

## Based on trophy
### Size
	kruskal.test(dat$Length..bp. ~ Habit, data = dat)
	kruskal.test(dat$Length..bp. ~ Habitat, data = dat)
	kruskal.test(dat$Length..bp. ~ Reproduction, data = dat)
	kruskal.test(dat$Length..bp. ~ Intermediate.host, data = dat)
### %GC
	kruskal.test(dat$Genome.GC. ~ Habit, data = dat)
	kruskal.test(dat$Genome.GC. ~ Habitat, data = dat)
	kruskal.test(dat$Genome.GC. ~ Reproduction, data = dat)
	kruskal.test(dat$Genome.GC. ~ Intermediate.host, data = dat)
### PCG Size
	kruskal.test(dat$PCG.Length ~ Habit, data = dat)
	kruskal.test(dat$PCG.Length ~ Habitat, data = dat)
	kruskal.test(dat$PCG.Length ~ Reproduction, data = dat)
	kruskal.test(dat$PCG.Length ~ Intermediate.host, data = dat)
### PCG % of Genome
	kruskal.test(dat$PCG...of.Genome ~ Habit, data = dat)
	kruskal.test(dat$PCG...of.Genome ~ Habitat, data = dat)
	kruskal.test(dat$PCG...of.Genome ~ Reproduction, data = dat)
	kruskal.test(dat$PCG...of.Genome ~ Intermediate.host, data = dat)
### Genome GC Skew
	kruskal.test(dat$Genome.GC.Skew ~ Habit, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Habitat, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Reproduction, data = dat)
	kruskal.test(dat$Genome.GC.Skew ~ Intermediate.host, data = dat)
### PCG GC Skew
	kruskal.test(dat$PCG.GC.Skew ~ Habit, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Habitat, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Reproduction, data = dat)
	kruskal.test(dat$PCG.GC.Skew ~ Intermediate.host, data = dat)
### Codon Position 3 Skew
	kruskal.test(dat$COX1.Codon.POS.3.GC.Skew ~ dat$Habit)
	kruskal.test(dat$COX1.Codon.POS.3.GC.Skew ~ dat$Habitat)
	kruskal.test(dat$COX1.Codon.POS.3.GC.Skew ~ dat$Reproduction)
	kruskal.test(dat$COX1.Codon.POS.3.GC.Skew ~ dat$Intermediate.host)

## Kruskal-Wallis post-hoc tests
	dunnTest(dat$PCG.GC.Skew ~ dat$Habit, data = dat)

## Additional testing
	kruskal.test(sub.dat1$Length..bp. ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$Length..bp. ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$Length..bp. ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$Length..bp. ~ Intermediate.host, data = sub.dat1)

	kruskal.test(sub.dat1$Genome.GC. ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC. ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC. ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC. ~ Intermediate.host, data = sub.dat1)

	kruskal.test(sub.dat1$PCG.Length ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.Length ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.Length ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.Length ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$PCG...of.Genome ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$PCG...of.Genome ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$PCG...of.Genome ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$PCG...of.Genome ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$Genome.GC.Skew ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC.Skew ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC.Skew ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC.Skew ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$PCG.GC.Skew ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.GC.Skew ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.GC.Skew ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.GC.Skew ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$avg.dn.rates ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$avg.dn.rates ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$avg.dn.rates ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$avg.dn.rates ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$avg.ds.rates ~ Habit, data = sub.dat1)
	kruskal.test(sub.dat1$avg.ds.rates ~ Habitat, data = sub.dat1)
	kruskal.test(sub.dat1$avg.ds.rates ~ Reproduction, data = sub.dat1)
	kruskal.test(sub.dat1$avg.ds.rates ~ Intermediate.host, data = sub.dat1)
	
	kruskal.test(sub.dat1$Length..bp. ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC. ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.Length ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$PCG...of.Genome ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$Genome.GC.Skew ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$PCG.GC.Skew ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$avg.dn.rates ~ Subclass, data = sub.dat1)
	kruskal.test(sub.dat1$avg.ds.rates ~ Subclass, data = sub.dat1)
