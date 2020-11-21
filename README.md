# ROCAPS
#### Supplementary data and scripts used in the article "Bayesian optimization of multivariate genomic prediction models based on secondary traits for improved accuracy gains and phenotyping costs" [[bioRxiv](biorxiv.org/content/10.1101/2020.11.20.390963v1)]
#### Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)
#### Date: 2020/11/13  (Last update: 2020/11/20)
Here, I explain the structure of this repository.

----------

- `ROCAPS`
	- `data.zip`: This folder contains datasets analyzed in this study. Please download and then decompress.
		- `extra`: This folder contain datasets other than genotype and phenotype.
			- `0.2_accession_group_information.csv`: cluster IDs and group IDs for each accession used in this study
		-  `genotype`: The marker genotype used in this study.
			-  `0.1_additive_genetic_relationship_matrix` (`RData`, `csv`): additive genomic relationship matrix used in this study
			-  `0.1_genotype_map.csv`: genetic map for marker genotype
			-  `0.1_genotype_matrix` (`RData`, `csv`): marker genotype data used in this study
		-  `phenotype`: The phenotypic data (real data & simulated data) in this study.
			-  `0.1_phenotype_selected` (`RData`, `csv`): real phenotypic values of selected traits used in this study
			-  `0.2_all_simulated_values_of_target_traits.RData`: an array of genotypic, residual, and phenotypic values of simulated target traits for 8 conditions (4 scenarios x 2 heritabilities)
			-  `0.2_observed_and_simulated_phenotype_for_ROCAPS` (`RData`, `csv`): all real target traits and simulated secondary traits used in this study
			-  `0.2_20_simulated_target_values_with2_herits_4_scenarios`: information on genotypic, residual, and phenotypic values of simulated target traits for 8 conditions (4 scenarios x 2 heritabilities) including variance-covariance structure
	- `raw_data.zip`: This folder contains raw datasets analyzed in this study. Please download and then decompress.
		- `RiceDiversity`: This folder contain original datasets downloaded from the website of "Rice Diversity" project (http://www.ricediversity.org/data/)
		-  `genotype`: The marker genotype used in this study.
			-  `geno.csv`: original marker genotype data downloaded from the website of "Rice Diversity" project
		-  `phenotype`: The phenotypic data used in this study.
			-  `pheno.csv`: original phenotypic data (real data) downloaded from the website of "Rice Diversity" project
	-  `scripts`: This folder contains scripts in the `R` language used in this study.
		-  `0.1_Rice_ROCAPS_data_modification.R`: Extract some traits used in this study from original phenotypic data, remove missing values, and calculate genomic relationship matrix
		-  `0.2_Rice_ROCAPS_generate_pseudo_data_of_plant_height_in_early_growth_stage.R`: Simulate growth curve using sigmoid curve, divide population into three groups, save all possible combinations of secondary traits, define cost, and simulate target traits 
		-  `1.0_Rice_ROCAPS_data_scaled_imputation.R`: function to scale and impute phenotypic data
		-  `1.1_Rice_ROCAPS_MTM_prediction.R`: function to perform multivariate GP using MTM
		-  `1.1.0_Rice_ROCAPS_MTM2.R`: function that modified MTM to control the message shown in the console
		-  `1.2_Rice_ROCAPS_rrBLUP_prediction.R`: function to perform GP by regarding secondary traits as fixed effects (not used in this study)
		-  `1.3_Rice_ROCAPS_RF_prediction.R`: function to perform GP by using Random Forest with explanatory variables of secondary traits and eigen vectors of marker genotype (not used in this study)
		-  `1.4_Rice_ROCAPS_Gaussian_process_for_BO.R`: function to perform Gaussian process (exploration step) for Bayesian optimization
		-  `1.5_Rice_ROCAPS_Acquisition_func.R`: function to compute acquisition function in Bayesian optimization
		-  `1.6_Rice_ROCAPS_convex_check.R`: function to judge whether the checked points are on the Pareto frontier or not
		-  `1.7_Rice_ROCAPS_pforeach`: function to enable us to perform parallelization with multiple cores
		-  `2.1_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation.R`: main script to optimize multivariate GP for gains in accuracy and phenotyping costs of secondary traits via Bayesian optimization using functions above
		-  `2.1.1_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation_for_terminal.R`, `2.1.2_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation_for_terminal_from_midstream.R`, `2.1.3_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Real_for_terminal.R`, `2.2.1_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation_for_terminal_run.txt`: scripts to run `2.1_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation.R` in different parameters with parallelization in terminal
			-  `3.1_Rice_ROCAPS_Summary_the_results_of_ROCAPS.R`: summary and plot the results obtained from `2.1_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Simulation.R`
