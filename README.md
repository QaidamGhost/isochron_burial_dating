**MATLAB scripts for calculating isochron burial age using paired cosmogenic nuclides.**

**June 10, 2022**

The MATLAB script "york.m" is ported from "York.R" which was written in R language and included in a toolbox called IsoplotR by Vermeesch, 2018, the "york_fixed_intercept.m" is modified after an unpublished R script written by Pieter Vermeesch (personal communication), and the "production_rate.m" is modified after a MATLAB script "CNP.m" developed by Chia-Yu Chen, Richard Ott, Erica Erlanger, Maarten Lupker and Yanyan Wang (Lupker et al., 2012). Other .m-files in this directory are written by Yizhou Yang. **See file "NOTICE" for copyright and licensing notice pertaining to all .m-files in this directory.**

## List of .m-files
 - isochron_burial_age: An iteration process to calculate isochron line for burial dating following Erlanger et al., 2012; Erlanger, 2010; Granger, 2014
 - KDE: Find the most probable value and 1 sigma absolute error for the observed simulations of an random variable after establishing the probability density function by a "nonparametric estimation" approach knows as kernel density estimation (KDE).
 - Npb_depth: Calculate the maximum estimation of the post-burial concentrations with 1 sigma error at the sampling location and penetration depth if the maximum exposure age (and additional erosion rate) of the surface is constrained.
 - plot_isochron: Plot isochron line, production ratio line, and errorbars of samples
 - production_rate: Calculate 10Be or 26Al Production Rates on the surface or at a given depth from Spallation and Muons.
 - remove_outliers: Statistically remove the reworked clast(s) following Odom, 2020
 - simple_burial_age: Calculate the simple burial age of each sample following the equation 22 in Granger, 2014.
 - york: A "york fit" script modified after york.R in IsoplotR by Vermeesch, 2018
 - york_fixed_intercept: A "york fit with a fixed intercept" script modified after an unpublished R script written by Pieter Vermeesch (personal communication)

## Instruction
 1. Define the variables, which are the arguments in "isochron_burial_dating.m", with your values.
 2. Create a file "expo_age.mat" (and an optional "e.mat") in the folder where the .m-files are. See the instruction of these two files in "Npb_depth.m" for the creation.
 3. Run "isochron_burial_age(data,init_Rinh,limit,source_lat,source_elv,measured_lat,measured_elv,z,rho,alpha);" in the command window.
 4. Wait for the plots and values of burial age in the command window.

## References
 - Balco, Greg, and Charles W. Rovey. "An isochron method for cosmogenic-nuclide dating of buried soils and sediments." American Journal of Science 308.10 (2008): 1083-1114.
 - Erlanger, Erica D. "Rock uplift, erosion, and tectonic uplift of South Africa determined with cosmogenic aluminum-26 and beryllium-10." Ph. D. Thesis (2010).
 - Erlanger, Erica D., Darryl E. Granger, and Ryan J. Gibbon. "Rock uplift rates in South Africa from isochron burial dating of fluvial and marine terraces." Geology 40.11 (2012): 1019-1022.
 - Gibbon, Ryan J., et al. "Early Acheulean technology in the Rietputs Formation, South Africa, dated with cosmogenic nuclides." Journal of Human Evolution 56.2 (2009): 152-160.
 - Granger, D. E. "Cosmogenic nuclide burial dating in archaeology and paleoanthropology." (2014): 81-97.
 - Lupker, Maarten, et al. "10Be-derived Himalayan denudation rates and sediment budgets in the Ganga basin." Earth and Planetary Science Letters 333 (2012): 146-156.
 - Mahon, Keith I. "The New “York” regression: Application of an improved statistical method to geochemistry." International Geology Review 38.4 (1996): 293-303.
 - Titterington, D. Michael, and Alex N. Halliday. "On the fitting of parallel isochrons and the method of maximum likelihood." Chemical Geology 26.3-4 (1979): 183-195.
 - Vermeesch, Pieter. "IsoplotR: A free and open toolbox for geochronology." Geoscience Frontiers 9.5 (2018): 1479-1493.
 - William III, E. Odom. Dating the Cenozoic incision history of the Tennessee and Shenandoah Rivers with cosmogenic nuclides and 40Ar/39Ar in manganese oxides. Diss. Purdue University Graduate School, 2020.
 - York, Derek. "Least squares fitting of a straight line with correlated errors." Earth and planetary science letters 5 (1968): 320-324.
 - York, Derek, et al. "Unified equations for the slope, intercept, and standard errors of the best straight line." American journal of physics 72.3 (2004): 367-375.
