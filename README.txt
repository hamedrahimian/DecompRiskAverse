Decomposition algorithms and scenario generation codes in C++

1. Subdirectories "D1C" to "D5S" contains the decomposition algorithms as introduced in the paper.
	1.a. The inputs for each decomposition algorithm code are files "input1_%d.txt", "input2.txt", "DemPortion.txt", "config.txt", 
		and a scenario file "scenFac_%d.txt". 
	     Variable No_scen should also be set to the number of children per node and must match %d in the scenario file name. 
	1.b. The outputs of each decomposition algorithm code are the log of the iterations and the CPU time.
2. ScenGen.cpp to generate scenario files 
	1.a. The inputs are variables numStage and numScen. By default, they are set to 4 and 25, respectively.  
	1.b. The output is "scenFac_%d.txt", where d is the number of children per node.
	



