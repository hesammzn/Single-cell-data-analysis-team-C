The files in this folder were executed in this order:

1- fastp_bash.sh:
	This script was used to trim the fastq files. The executable file for fastp (version 0.24.0) is needed to run this script.

2- salmon_index.sh:
	This code was used to generate gene indices from transcript.fa file. The executable file for salmon (version 1.10.0) and files inside the library folder (comes with the salmon compressed file) are needed to be at the same directory as the terminal and .fa file for this code to run.

3- salmon_map.sh:
	This script was used to map the reads to index built in the step 2.

	



Other files in this folder:

accessions:
	Cell accession numbers.
