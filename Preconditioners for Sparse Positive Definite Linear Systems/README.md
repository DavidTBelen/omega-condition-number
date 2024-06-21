# Organization of the folder:

1. run_file.m 

	- Runs the main experiment.
	- It needs to specify the matrices data sets.
	- Saves an output file with the results of the experiment.

2. performance_profile_file.m

	- Plots the performance profile in Section 4.1 in the paper by loading the results file.
	- It calls for "tps_performance.m" and "perf.m"

3. table_file.m
	
	- Generates the tables from Section 4.1 in the paper by loading the results file.
	- Choose type for obtaining the desired table (time, iterations, residuals,...)

4. selected_matrices.m
	
	- Takes a subset of test matrices from the USF repository
	  and removes the ones that lead to "trivial problems"
	- Output is a file "selected_matrices_list.mat" that can be used as data for the "run_file"

5. folder "data" contains the stored data from the experiments
	 
	- Folder "suitesparse_matrices" contains .mat files of matrices.

6. folder "Auxiliary_files" contains different functions:

	- Folder "Preconditioners" contain the codes for computing the diagonal and diagonal + block triangular preconditioner
