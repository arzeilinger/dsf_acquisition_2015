Analysis of data from DSF (diffusible signalling factor) acquisition experiment with Xylella fastidiosa.

R_analyses/hierarchical_transmission_model.R includes code for Bayesian hierarchical model for analysis of transmission data using nimble package.

R_analyses/DSF_acquisition_analysis.R includes code using a series of GLMM models.

<b>qPCR Calculations</b>

R_functions/qpcrFunctions.R provides functions for merging Cq/N0 values from LinRegPCR with sample names and serial dilution data. Xylella CFUs can then be calculated from the standard curve and N0 estimates.

R_analyses/qpcr_calculates.R uses calculates Xylella CFU from qPCR output for the DSF acquisition experiment.
