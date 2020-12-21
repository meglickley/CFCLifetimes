Lickley et al. (in Revisions) Joint inference of CFC lifetimes and banks suggests previously unidentified emissions

Authors: Megan Lickley, Sarah Fletcher, Matt Rigby and Susan Solomon
Corresponding Author: Megan Lickley; mlickley@mit.edu


Code for generating figures in main text. 

Contained in Code Directory

submodelRFandDE.m: matlab script that simulates the bottom up accounting of banks by equipment type and generates joint prior distributions of RF (release fraction) and DE (direct emissions)

MainScript_JointInference.m: matlab script that runs the BPE model. 

simulation_model.m: matlab function that is called from MainScript_multi.m

Bank_Emiss.m:  matlab function that is called from simulation_model.m

MakeFigures.m: matlab script that makes the figures after the prior and posterior samples have been generated for each scenario. 

Figure_Data: a directory containing data from previous published assessments, the WACCM model, and results from simulations of the BPE under various assumptions that are used to create figures in the main text  

boundedline.m: a function used to create the figures. 

Input: a directory containing input data divided by CFC11, CFC12 and CFC113, as well as SPARC_lifetimes. This directory will also contain output from submodelRFandDE.m and from estimating_prod.m.  Each of the directories within the Input directory for each gas contains the following: 

- a5_na5.mat: a file containing production data used in MainScript_multi.m
- AFEAS_cfcXXproduction.mat: a file containing AFEAS production data by equipment type
- WMO2002.mat: a file containing production data used in MainScript_multi.m
- wmo2018.mat: a file containing atmospheric concentrations 



%%%%%%%%%%%%%%%%%%  INSTRUCTION FOR RUNNING THE BPE MODEL %%%%%%%%%%%%%%%%%%%%%%

Step 1:  Open submodelRFandDE.m script.  Change HomeDir to specify location of file. Run this script.  This will generate joint priors for RF and DE that is then stored in the Input Folder under the corresponding gas' directory. This should take approximately 20 minutes on a “normal” desktop computer. 

Step 2:  Open estimating_prod.m script.  Change HomeDir to specify location of file. Run this script.  This will generate prior samples of production.  This should take approximately 10 minutes. 

Step 3:  Open MainScript_JointInference.m.  Change HomeDir to specify location of file.  Run the following scenarios by changing the variables under‘Scenarios’ on lines 36 and 37
	- unreportedProd = 1, jointInference = 1; 
	- unreportedProd = 0, jointInference = 1; 	
	- unreportedProd = 1, jointInference = 0; 
Each simulation should take 20-30 minutes on a normal desktop computer. If the model is running correctly then figures should look approximately like the example output figures in the Output directory:  BPE_modeloutput_Fig1.jpg, and BPE_modeloutput_Fig2.jpg


Step 4:  Open MakeFigures.m and change HomeDir to specify location of the Code folder.  Run the script to generate figures and populate table values.  All figures will be saved in the FigureFolder directory. Table values will be printed in the Command Window. Creating each figure should take approximately 1 minute on a “normal” desktop computer. 