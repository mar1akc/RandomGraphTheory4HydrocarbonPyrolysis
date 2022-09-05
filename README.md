# Random Graph Theory for Hydrocarbon Pyrolysis
The configuration model and generating function formalism proposed by Newman, Strogats and Watts (2001) 
(https://journals.aps.org/pre/pdf/10.1103/PhysRevE.64.026118)
can be used to predict molecule size distribution in hydrocarbon pyrolysis 
(Dufour-Decieux, Moakler, Cameron, Reed, 2022, https://arxiv.org/abs/2205.13664).

This repository contains a collection of Matlab codes allowing you to reproduce Figures 4 and 5 in the main text and S6, S7, and S8 in Supplementary Materials.

Figures 4, S6, and S7 display size distributions for small molecules 

(1) predicted by random graph theory from the degree distribution extracted from MD simulations (RGT), 

(2) predicted by random graph theory from the degree distribution predicted by the ten-reaction model (10RM+RGT),

(3) and extracted from MD simulations,

and Wasserstein W1 dstances between distributions [(1) and (3)] and [(2) and (3)].

To reproduce small molecule size distributions in Figures 4, S6, and S7, run SmallMolSizeDistr.m. 

To reproduce plots with Wasserstein distances, run W1dist4SmallMolSizeDistr.m
