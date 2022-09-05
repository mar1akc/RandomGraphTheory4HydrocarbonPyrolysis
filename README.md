# Random Graph Theory for Hydrocarbon Pyrolysis
The configuration model and generating function formalism proposed by Newman, Strogats and Watts (2001) 
(https://journals.aps.org/pre/pdf/10.1103/PhysRevE.64.026118)
can be used to predict molecule size distribution in hydrocarbon pyrolysis 
(*Dufour-Decieux, Moakler, Cameron, Reed*, 2022, https://arxiv.org/abs/2205.13664).

This repository contains a collection of Matlab codes allowing you to reproduce Figures 4 and 5 in the main text and S6, S7, and S8 in Supplementary Materials.

Figures 4, S6, and S7 display size distributions for small molecules <br>
(1) predicted by random graph theory from the degree distribution extracted from MD simulations (RGT), <br>
(2) predicted by random graph theory from the degree distribution predicted by the ten-reaction model (10RM+RGT),<br>
(3) and extracted from MD simulations (MD),<br>
and Wasserstein W1 dstances between distributions [(1) and (3)] and [(2) and (3)].

To reproduce small molecule size distributions in Figures 4, S6, and S7, run ***SmallMolSizeDistr.m***. <br>
Input data for ***SmallMolSizeDistr.m*** are found in folder ***Data***. <br>  
Functions ***H0distributions.m***, ***CompSizeDistr_pi.m*** and ***CompSizeDistr_derivatives.m*** are called by ***SmallMolSizeDistr.m***. <br>
Function ***H0distribution.m*** computes the probability distribution $(P_s)$ where $P_s$ is the probability for a randomly picked vertex to belong to a connected component of size $s$. <br>
Function ***CompSizeDistr_pi.m*** recasts the distribution $(P_s)$ to the distribution $(\pi_s)$ where $\pi_s$ is the probability for a randomly picked connected component to contain $s$ vertices. <br>
Function ***CompSizeDistr_derivatives.m*** is used to compute error bars for distribution (1). <br> 
Function ***SmallMolSizeDistr.m*** generates an input file ***W1data.mat*** for ***W1dist4SmallMolSizeDistr.m***. 

To reproduce plots with Wasserstein distances, run ***W1dist4SmallMolSizeDistr.m***
