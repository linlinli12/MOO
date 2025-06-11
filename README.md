# Digital Cousins: Simultaneous Optimization of One Model for BMP Signaling in Distant Relatives Reveals Essential Core

This repository contains the source code, simulation scripts, and optimization workflows supporting the publication:
Linlin Li^1 and Thembi Mdluli^1,  Gregery Buzzard^2, David Umulis^(1,3)
1.	Weldon School of Biomedical Engineering Purdue University
West Lafayette, ID 47906 
2.	Department of Mathematics, Purdue University West Lafayette, ID 47906
3.	Weldon School of Biomedical Engineering Purdue University
West Lafayette, ID 47906, dumulis@purdue.edu


Abstract:
Spatially distributed, nonuniform morphogen gradients play a crucial role in tissue organization during development across the animal kingdom. The Bone Morphogenetic Protein (BMP) pathway, a well-characterized morphogen system involved in dorsal-ventral (D-V) axis patterning, is conserved in both form and function across species such as Drosophila and zebrafish. This study investigates whether a shared core mathematical model—with a fixed network topology and constrained biophysical parameters—can explain BMP gradient formation in both species using species-specific gene expression profiles alone.

Using multi-objective optimization, we simultaneously fit a single model to BMP data from Drosophila and zebrafish. This enabled us to uncover essential features of the pathway that are conserved and highlight parameters responsible for interspecies differences. Our findings show that only two parametric changes are sufficient for the model to accurately reproduce experimental BMP gradients in both species. Rather than seeking predictive generalization, this cross-species modeling framework aims to enhance the interpretability and biological relevance of computational models of conserved developmental pathways.

## Code Availability
All MATLAB code used in this study is available at [GitHub link or repository DOI]. The repository contains scripts for parameter screening, front analysis, and result evaluation relevant to the computational modeling described in the manuscript.

The folder structure includes:

ExpData/: Experimental data used for model fitting and validation.

Function/: Utility functions used across analyses.

Main scripts:

Front_analysis_WT.m and Front_analysis_chordinlike.m: Simulations and analyses of wild-type and chordin-like scenarios.

Front_analysis_find_best_parameter.m: Identifies optimal parameters based on front position metrics.

Parameter_screen_DegreeFreedom.m: Main code for running the parameter screen

Parameter_screen_localsensitivity_updated.m: Computes local parameter analysis around best parameter set.

parameter_gen.m: Generates parameter sets using Latin Hypercube sampling.

result_analysis_Linlin_multiFolder_degree_Chorinlike.m and result_analysis_Linlin_multiFolder_degree_final.m:  Aggregate and visualize results across multiple simulation for chordin-like and WT results, respectively.
best_parm_200x5000.mat: Stored best-fitted parameter sets from the optimization screen.
