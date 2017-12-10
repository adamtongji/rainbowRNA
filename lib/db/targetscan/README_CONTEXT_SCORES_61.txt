INTRODUCTION

The Perl script targetscan_60_context_scores.pl calculates context+ scores for a set of miRNA targets predicted by targetscan_60.pl. 
The TargetScan 6.1 context+ score code produces essentially the same output as displayed at TargetScan.org Release 6.1.  
The only change in the calculation of context+ scores from Release 6.0 ro Release 6.1 is the inclusion of TA and SPS parameters 
for all miRNA families (rather than just highly conserved families).  Thus TA_SPS_by_seed_region.txt has been updated 
(along with targetscan_60_context_scores_output.txt, containing sample output) but all other files are the same.

The script takes four input files
	1) a miRNA file: a tab-delimited text file of mature miRNA information.  This is different from the file required by targetscan_60.pl.
	2) a UTR file: a tab-delimited multiple sequence alignment of the 3' UTRs of genes from the desired species
		which is the same as the input file for targetscan_60.pl.
	3) a PredictedTargets file: output from targetscan_60.pl.
	4) TA_SPS_by_seed_region.txt: Contains input parameters and must be in same directory as targetscan_60_context_scores.pl.

In this directory we are providing samples of all of the above files.

The sample files are: miR_for_context_scores.txt, UTR_sequences_sample.txt and targetscan_60_output.txt
If you wish to generate these files from the complete data (linked from the table on the Data Download page), 
	run the commands shown below to convert them to the correct format.


FILE FORMATS

The format of the input files is important for the script to work correctly. 

Each line of the miRNA mature sequence file (ex: miR_for_context_scores.txt) consists of 4 tab separated entries:
1) miRNA family ID: Name of the miRNA family
2) Species ID of this miRNA family (which should match species IDs in the UTR and predicted targets input files)
3) MiRBase ID: name of a mature miRNA sequence
4) Mature sequence: sequence of mature miRNA

To generate a file in this format from the complete "miR Family" file (miR_Family_Info.txt from the table 
on the Data Download page), run this command:
cut -f1,3,4,5 miR_Family_Info.txt > miR_for_context_scores.txt

Each line of the UTR alignment file (ex: UTR_sequences_sample.txt) consists of 3 tab separated entries
1) Gene symbol or transcript ID
2) Species ID (which should match species IDs in miRNA input file) 
3) Sequence 

To generate a file in this format from the complete "UTR Sequences" file (UTR_Sequences.txt from the table 
on the Data Download page), run this command:
cut -f1,4,5 UTR_Sequences.txt > UTR_sequences_sample.txt

Each line of the predicted targets file (ex: targetscan_60_output.txt) consists of 13 tab separated entries
(although not all fields are required)
1)  GeneID - name/ID of gene (from UTR input file)
2)  miRNA family_ID - name/ID of miRNA family (from miRNA input file)
3)  species ID - name/ID of species (from UTR input file)
4)  MSA start - starting position of site in aligned UTR (counting gaps) 
5)  MSA end - ending position of site in aligned UTR (counting gaps) 
6)  UTR start - starting position of site in UTR (not counting gaps) 
7)  UTR end - ending position of site in UTR (not counting gaps) 
8)  Group ID - ID (number) of site(s) (same gene, same miRNA) that overlap 
9)  Site type - type of site in this species (1a [7mer-1a; type 1], m8 [7mer-m8; type 2], or 8mer [type 3])
10) miRNA in this species - if "x", then this miRNA has been annotated in this species
11) Group type - type of this group of sites; if 'Site_type' in a 'Group_ID' is heterogeneous, "weakest" type of the group is used
12) Species in this group - list of species names/IDs in which this site is found
13) Species in this group with this site type - for hetergeneous groups only


EXECUTION

The script can be executed in 3 different ways:
1) Running the script without any arguments (./targetscan_60_context_scores.pl) will print out a help screen.
1) Running the script without the '-h' flag (./targetscan_60_context_scores.pl -h) will print out a formats of input files.
2) Running the script with input filenames and output file will perform the analysis. Ex:
	./targetscan_60_context_scores.pl miR_for_context_scores.txt UTR_Sequences_sample.txt targetscan_60_output.txt targetscan_60_context_scores_output.txt

OUTPUT FILES

In this folder is a sample output file called "targetscan_60_context_scores_output.txt".
The output file also contain several tab separated entries per line.

The sample output file has a headers that names each column:
1)  Gene ID - name/ID of gene (from UTR input file)
2)  Species ID - name/ID of species (from UTR input file)
3)  Mirbase ID - name of a mature miRNA sequence
4)  Site Type - type of site in this species (1a [7mer-1a; type 1], m8 [7mer-m8; type 2], or 8mer [type 3])
5)  UTR start - starting position of site in UTR (not counting gaps) 
6)  UTR end - ending position of site in UTR (not counting gaps) 
7)  3' pairing contribution - one component of context+ score (see Grimson et al., 2007 for details)
8)  local AU contribution - one component of context+ score (see Grimson et al., 2007 for details)
9)  position contribution - one component of context+ score (see Grimson et al., 2007 for details)
10) TA contribution - one component of context+ score (see Garcia et al., 2011 for details)
11) SPS contribution - one component of context+ score (see Garcia et al., 2011 for details)
12) context+ score - context+ score, the sum of the above five contributions and site-type contribution
	(-0.310 for 8mer, -0.161 for 7mer-m8, -0.099 for 7mer-1A)
13) context+ score percentile rank - percentage of sites for this miRNA with a less favorable context score
14) UTR region - subsequence of UTR used to show predicted consequential pairing
15) UTR-miRNA pairing - predicted consequential pairing of target region and miRNA; complementary bases indicated by bars
16) mature miRNA sequence - mature sequence of this Mirbase ID
17) miRNA family - name/ID of miRNA family
18) Group # - ID (number) of site(s) (same gene, same miRNA) that overlap 

NOTES

This script was designed on a Linux platform. While running this script on Windows or Mac platforms, make sure to call the native perl binary. 
This can be done explicitly by executing it as 'perl targetscan_60_context_scores.pl' or changing the first line of the script to point to the native binary. 


QUESTIONS/SUGGESTIONS:

Please direct all correpondence to wibr-bioinformatics@wi.mit.edu




TA contribution
SPS contribution
context+ score
context+ score percentile
UTR region
UTR-miRNA pairing
mature miRNA sequence
miRNA family
Group #
