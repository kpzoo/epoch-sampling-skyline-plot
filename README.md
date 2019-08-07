# epoch-sampling-skyline-plot
Code for the maximum likelihood solution to the epoch sampling skyline plot (ESP) method of Parag, du Plessis and Pybus

- the ESP is implemented in Matlab and processes data generated from R scripts
- utilises the phylodyn R package of Karcher et al and the shadedplot Matlab package

To use:
1) Generate coalescent tree data using one of the R scripts in the using phylodyn folder
2) Resulting csvs can be loaded into a folder of choice
3) Use a main Matlab function (e.g. batchSwitchBetaSimple) to run the ESP on the csvs
4) Outputs statistics and plots of estimated effective population sizes and sampling intensities



References:

ESP:

K. Parag, L. du Plessis, and O. Pybus. An integrated framework for the joint inference of demographic history and sampling intensity from genealogies or genetic sequences. BioRxiv, 686378, 2019.



Dependencies: 

phylodyn: M. Karcher, J. Palacios, S. Lan, et al. PHYLODYN: an R package for Phylodynamic Simulation and Inference. Mol. Ecol. Res, 17: 96–100, 2017.

shadedplot: D. Van Tol. Shaded area plot (https://uk.mathworks.com/matlabcentral/fileexchange/18738-shaded-area-plot), MATLAB Central File Exchange, 2005.
