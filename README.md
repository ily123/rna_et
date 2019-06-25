This code repository accompanies the manuscript

	*An Evolutionary Trace method defines functionally
	important bases and sites common to RNA families*

There are two components of interest

1. ET program (ET_RNA_ms/et_bin/wetc)
	stand-alone C binary that produces ET ranking given a MSF alignment
			
	To run:
	1) In terminal, navigate to ET_RNA_ms/et_bin/
	2) Ensure the program runs by typing "./wetc - help"
	3) To score an RNA aligment, type
	
		./wetc -p [aln_path] -x [query_seq_name] -dnaet -realval -o [output_name]

		-p 		path to MSF alignment
		-x 		name of the query sequence as it appears in the alignment
		-dnaet		flag to enable non-protein scoring
		-realval	flag to run rvET scoring
		-o		base name for output files (can pipe it directly
          into dir '-o dir_path/basename')

2. Example of RNAse P analysis, wrapped in MATLAB.

   To run:
   1. Launch MATLAB
   2. Navigate MATLAB console to ET_RNA_ms/rnasep_example/
   3. Execute script 'rnasep_demo.m' by typing 'rnasep_demo'
   4. The script will:

      a. Load msf alignment of RNAse P homologues

      b. Load RNAse P structure 3q1q (chain B)

      c. Calculate adjacency matrix for 3q1q (chain B)

      d. Profile sequence of chain B to the MSF alignment as query 

      e. Trace the profile alignment
           (outputs are saved to ../et_temp/, then mv'ed
            to ./results) 

      f. Calculate random expected clustering for adj matrix

      g. Calcuate clustering by ET nucleotides (across bins 0-100).

      h. Calculate overlap of ET nucleotides with functional sites

      i. Print summary tables and figures to ./results and ./visualization

   5. The necessary inputs for this script are in ./input_data/
      and include:

      1. MSF alignment of RNAse P homologues
      2. 3Q1Q pdb file

3. Possible Errors:
    1. Make sure gaps are denoted as '.' in the MSF alignment
    2. If using wetc standalone, make sure the file path in -o
       is not too long (will throw error)


Bonus:

  Ribosome Example (ET_RNA_ms/ribosome_23S_example)
  example of ET application using the 23S rRNA alignment,
  including the ET trace, and z-score computations
  (Several intermediate steps here are already pre-computed)

  To run:
  1. Launch matlab
  2. Navigate in matlab console to ET_RNA_ms/ribosome_23S_example
  3. Type 'ls' (no quotes) to examine contents of folder, there are 3 scripts of interest
    trace_alignment.m     traces the 16S alignment
    measure_et_structural_clustering  calculates ET clustering in the 16S structure
    measure_et_prediction_accuracy    calculates ET overlap with know functional sites in 16S
  4. Execute these 3 scripts in order
  5. Results and plots will be in the results/ and visualizaion/ folders

Dependencies: The code was tested on MATLAB v2017b with Bioinformatics toolbox.
