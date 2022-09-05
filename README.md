# Random Graph Theory for Hydrocarbon Pyrolysis
The configuration model and generating function formalism proposed by Newman, Strogats and Watts (2001) 
(https://journals.aps.org/pre/pdf/10.1103/PhysRevE.64.026118)
can be used to predict molecule size distribution in hydrocarbon pyrolysis 
[*Dufour-Decieux, Moakler, Cameron, Reed*, 2022, https://arxiv.org/abs/2205.13664].

This repository contains a collection of Matlab codes allowing you to reproduce Figures 4 and 5 in the main text and S6, S7, and S8 in Supplementary Materials.

Figures 4, S6, and S7 display **size distributions for small molecules** <br>
(1) predicted by random graph theory from the degree distribution extracted from MD simulations (RGT), <br>
(2) predicted by random graph theory from the degree distribution predicted by the ten-reaction model (10RM+RGT),<br>
(3) and extracted from MD simulations (MD data),<br>
and Wasserstein W1 dstances between distributions [(1) and (3)] and [(2) and (3)].

* To reproduce small molecule size distributions in Figures 4, S6, and S7, run ***SmallMolSizeDistr.m***. <br>
Input data for ***SmallMolSizeDistr.m*** are found in folder ***Data***. There are two input files: ***Degrees_predictions_10reac.csv*** and ***DegreeAndMolSizeMDdata.mat***. The total number of datasets available is 17, the first 14 of which were used for the main study conducted in [*Dufour-Decieux, Moakler, Cameron, Reed*, 2022, https://arxiv.org/abs/2205.13664*] while the last three with 320 carbons each were generated for the study of the size effect.<br>
The file ***Degrees_prediction_10reac.csv*** contains degree counts predicted by the ten-reaction model, as well as  the numbers of carbons, temperatures, and H/C ratio data. In order to convert the degree counts to degree distributions, we divide these counts by the numbers of carbons in corresponsing systems also given in this file. ***Degrees_predictions_10reac.csv*** is converted to data arrays by function ***read_data(filename)***. <br>
The file ***DegreeAndMolSizeMDdata.mat*** is a structure consisting of two 3-by-17 cell arrays: ***DegreeDistribution_cell*** and ***MolSizeDistribution_cell***:<br>
d = load(fname_DegAndMolSize);<br>
dd = d.DegreeDistribution_cell;<br>
% dd = 3-by-17 cell array<br>
% dd{1,j} = string with the name of the dataset j<br>
% dd{2,j} = degree distribution for dataset j extracted from MD simulations<br>
% dd{3,j} = standard deviations for degree distribution j extracted from MD simulations<br>
sd = d.MolSizeDistribution_cell; <br>
% sd = 3-by-17 cell array <br>
% sd{1,j} = string with the name of the dataset j <br>
% sd{2,j} = molecule size distribution for dataset j extracted from MD simulations <br>
% sd{3,j} = standard deviations for molecule counts for dataset j  <br>
Functions ***H0distributions.m***, ***CompSizeDistr_pi.m*** and ***CompSizeDistr_derivatives.m*** are called by ***SmallMolSizeDistr.m***. <br>
Function ***H0distribution.m*** computes the probability distribution $(P_s)$ where $P_s$ is the probability for a randomly picked vertex to belong to a connected component of size $s$. <br>
Function ***CompSizeDistr_pi.m*** recasts the distribution $(P_s)$ to the distribution $(\pi_s)$ where $\pi_s$ is the probability for a randomly picked connected component to contain $s$ vertices. <br>
Function ***CompSizeDistr_derivatives.m*** is used to compute error bars for distribution (1). <br> 
Function ***SmallMolSizeDistr.m*** generates an input file ***W1data.mat*** for ***W1dist4SmallMolSizeDistr.m***. 

* To reproduce plots with Wasserstein distances, run ***W1dist4SmallMolSizeDistr.m***

Figures 5 and S8 in [*Dufour-Decieux, Moakler, Cameron, Reed*, 2022, https://arxiv.org/abs/2205.13664*] display the **largest molecule size distributions**<br> 
(4) extracted from MD data (MD data) and ,br>
(5) obtained using random graph sampling for degree distributions generated using the ten-reaction model (10RM+RGS). <br>
To reproduce the histograms depicting the largest molecule size distributions run the following sequence of codes:<br>
* Run ***SampleConfigurationModel.m***. It generates Nsamples=10000 samples of random graphs for each degree distribution contained in the input file ***Data/Degrees_predictions_10reac.csv*** and creates output file ***Data/largest_mol_hist_data_10RM_RGS.mat***. File ***Data/largest_mol_hist_data_10RM_RGS.mat*** contains 3-by-17 cell array of the following format: <br>
data_10RM_RGT = load("Data/largest_mol_hist_data_10RM_RGS.mat");
d2 = data_10RM_RGT.largest_mol_hist_data_10RM_RGS; <br>
% d2{1,j} = string with the name of the dataset j<br>
% d2{2,j} = array of counts of random graph samples with largest connected component of size s<br>
% d2{3,j} = array of values s spanning the range of sizes of largest connected components<br>
* Run ***LargestMolHistograms.m*** to generate the histograms and compute the means and the standard deviations for the largest molecule size distributions. Input files: ***Data/largest_mol_hist_data_10RM_RGS.mat*** and ***Data/largest_mol_hist_data_MD.mat***. File ***Data/largest_mol_hist_data_MD.mat*** contains a 3-by-17 cell array of the same format as ***Data/largest_mol_hist_data_10RM_RGS.mat*** except for the largest molecula size distributions in it are extracted from MD data.
