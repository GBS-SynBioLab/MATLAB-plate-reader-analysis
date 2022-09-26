# MATLAB-plate-reader-analysis
Script to analyse plate-reader data (OD and fluorescence) in MATLAB

## Instructions:
- Data from "AB20200818_BW_409_410.xlsx" have been previously exported into MATLAB data "AB20200818_BW_409_410.mat" consisting of 2 matrices: OD700 and GFP
- Run "AB20200818_BW_409_410.m" in MATLAB to analyse the data from Excel file "AB20200818_BW_409_410.xlsx", the "barweb.m" and "stdshade.m" need to be present in the same folder or in an active MATLAB path.

## Script structure:
1) Define bacterial culture condition indexes in triplicate (wells that they correspond to in the 96-well plate)
2) Import data from "AB20200818_BW_409_410.mat" file
3) Subtract background OD from the M9 medium
4) Subract background fluorescence from the WT cells
5) Smooth the OD and GFP data with the Smoothing Spline MATLAB function
6) Calculate the growth rate and GFP production rate per cell according to Ceroni et al. 2015^{*}
7) Plot bar graph of growth rate and GFP at the timepoint of interest
8) Plot OD and growth rate over time
9) Plot GFP, GFP per cell and GFP production rate per cell over time
10) Plot the static input-output curves of the GFP vs inducer concentration
11) Plot bar graph of the maximum growth rate of each cell condition


^{*} Ceroni, Francesca, et al. "Quantifying cellular capacity identifies gene expression designs with reduced burden." Nature methods 12.5 (2015): 415-418.
