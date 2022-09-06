**MATLAB scripts for calculating isochron burial age using paired cosmogenic nuclides.**

**June 10, 2022**

The MATLAB script "york.m" is ported from "York.R" which was written in R language and included in a toolbox called IsoplotR by Vermeesch, 2018, the "york_fixed_intercept.m" is modified after an unpublished R script written by Pieter Vermeesch (personal communication), and the "production_rate.m" is modified after a MATLAB script "CNP.m" developed by Maarten Lupker, Chia-Yu Chen, Richard Ott, Erica Erlanger, and Yanyan Wang (Lupker et al., 2012). Other .m-files in this directory are written by Yizhou Yang. **See file "NOTICE" for copyright and licensing notice pertaining to all .mat and .m-files in this directory.**

## List of .mat and .m-files
 - consts.mat: constants for calculation, including Stone's scaling parameters, mean life of 10Be and 26Al, step for KDE, simalation times, SLHL production rates, limit for iterations, and attenuation lengths of nuclide and muons.
 - custom_mat_files: Generate custom "e.mat" and "expo_age.mat" before maximum estimation if erosion rates and exposure age have already been known.
 - isochron_burial_age: An iteration process to calculate isochron line for burial dating following Erlanger et al., 2012; Erlanger, 2010; Granger, 2014
 - KDE: Find the most probable value and 1 sigma absolute error for the observed simulations of an random variable after establishing the probability density function by a "nonparametric estimation" approach knows as kernel density estimation (KDE).
 - Npb_depth: Calculate the maximum estimation of the post-burial concentrations with 1 sigma error at the sampling location and penetration depth if the maximum exposure age (and additional erosion rate) of the surface is constrained.
 - plot_isochron: Plot isochron line, production ratio line, and errorbars of samples
 - production_rate: Calculate 10Be or 26Al Production Rates on the surface or at a given depth from Spallation and Muons.
 - remove_outliers: Statistically remove the reworked clast(s) following Odom, 2020
 - simple_burial_age: Calculate the simple burial age of each sample.
 - york: A "york fit" script modified after york.R in IsoplotR by Vermeesch, 2018
 - york_fixed_intercept: A "york fit with a fixed intercept" script modified after an unpublished R script written by Pieter Vermeesch (personal communication)

## Instruction
### Preparation
1.  The scripts assume the ratio of surface production rates as  **6.61**  which is indicated by the adopted mean life of 10Be and 26Al.  **Try to modify the production rates and mean life in the "consts.mat" and/or the codes in the "production_rate.m" if you want to switch to a different value of that ratio, or the results will show some inaccuracy in burial ages.**
2. Make sure that all [DEPENDENCIES](https://github.com/QaidamGhost/isochron_burial_dating#software-version-and-dependency) are all promised.
3. Create the file "expo_age.mat" (and an optional "e.mat") in the root directory. See the instruction of these two files in "Npb_depth.m" and "custom_mat_files" for the creation.
4. Check the values restore in the "consts.mat". If you want to use some other values, just replace them in the .mat-file before running the scripts.
### Calculation
1.  Define the variables, which are the arguments in "isochron_burial_dating.m", with your values. See "arguments" in "isochron_burial_dating.m" for further instructions.
```
    data.x=[XXXXXX,XXXXXX,...,XXXXXX];
    data.dx=[XXXXXX,XXXXXX,...,XXXXXX];
    data.y=[XXXXXX,XXXXXX,...,XXXXXX];
    data.dy=[XXXXXX,XXXXXX,...,XXXXXX];
    shielding_factor=XXXXXX;
    init_Rinh=XXXXXX;
    source_lat=XXXXXX;
    source_elv=XXXXXX;
    measured_lat=XXXXXX;
    measured_elv=XXXXXX;
    z=XXXXXX;
    rho=XXXXXX;
    alpha=XXXXXX;
    option.flag2=1; % optional
```
2.  Run "isochron_burial_age.m" in the command window.
 ```
	isochron_burial_age(data,init_Rinh,source_lat,source_elv,measured_lat,measured_elv,shielding_factor,z,rho,alpha,option);
```
### Output
1. Print the simple burial ages and their uncertainties of each sample.
2. Print the isochron burial ages, MSWD, and its uncertainties and plot the isochron line with measured data, linearized data, and reworked clasts.
3. If the final intercept in step 2 is less than zero, the scripts will give minimum and maximum estimations of the burial age. 
3.1. Print the minimum isochron burial ages, MSWD, and its uncertainties and plot the isochron line with measured data, linearized data, and reworked clasts.
3.2. Print the maximum isochron burial ages, MSWD, and its uncertainties and plot the isochron line with measured data, linearized data, and reworked clasts.

## Software Version and Dependency
 - Matlab 2021b or higher version (Recommend).
 - Optimization Toolbox 9.2 or higher version.
 - Statistics or Machine Learning Toolbox 12.2 or higher version.
 - Curve Fitting Toolbox 3.6 or higher version.

## Platforms for Tests
 - Intel i9-9900K / HP OMEN by HP Obelisk Desktop 875-1xxx / Windows 10 21H2
 - AMD RYZEN R5-3600x / MSI B450m MORTAR / Windows 10 21H2

## References
 - Balco, Greg, and Charles W. Rovey. "An isochron method for cosmogenic-nuclide dating of buried soils and sediments." American Journal of Science 308.10 (2008): 1083-1114.
 - Braucher, R., et al. "Production of cosmogenic radionuclides at great depth: A multi element approach." Earth and Planetary Science Letters 309.1-2 (2011): 1-9.
 - Chmeleff, Jérôme, et al. "Determination of the 10Be half-life by multicollector ICP-MS and liquid scintillation counting." Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms 268.2 (2010): 192-199.
 - Erlanger, Erica D. "Rock uplift, erosion, and tectonic uplift of South Africa determined with cosmogenic aluminum-26 and beryllium-10." Ph. D. Thesis (2010).
 - Erlanger, Erica D., Darryl E. Granger, and Ryan J. Gibbon. "Rock uplift rates in South Africa from isochron burial dating of fluvial and marine terraces." Geology 40.11 (2012): 1019-1022.
 - Gibbon, Ryan J., et al. "Early Acheulean technology in the Rietputs Formation, South Africa, dated with cosmogenic nuclides." Journal of Human Evolution 56.2 (2009): 152-160.
 - Granger, Darryl E., and Paul F. Muzikar. "Dating sediment burial with in situ-produced cosmogenic nuclides: theory, techniques, and limitations." Earth and Planetary Science Letters 188.1-2 (2001): 269-281.
 - Granger, D. E. "Cosmogenic nuclide burial dating in archaeology and paleoanthropology." (2014): 81-97.
 - Hidy, Alan J., et al. "A geologically constrained Monte Carlo approach to modeling exposure ages from profiles of cosmogenic nuclides: An example from Lees Ferry, Arizona." Geochemistry, Geophysics, Geosystems 11.9 (2010).
 - Korschinek, Gunther, et al. "A new value for the half-life of 10Be by heavy-ion elastic recoil detection and liquid scintillation counting." Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms 268.2 (2010): 187-191.
 - Lal, Devendra. "Cosmic ray labeling of erosion surfaces: in situ nuclide production rates and erosion models." Earth and Planetary Science Letters 104.2-4 (1991): 424-439.
 - Lupker, Maarten, et al. "10Be-derived Himalayan denudation rates and sediment budgets in the Ganga basin." Earth and Planetary Science Letters 333 (2012): 146-156.
 - Mahon, Keith I. "The New “York” regression: Application of an improved statistical method to geochemistry." International Geology Review 38.4 (1996): 293-303.
 - Samworth, E. A., E. K. Warburton, and G. A. P. Engelbertink. "Beta Decay of the Al 26 Ground State." Physical Review C 5.1 (1972): 138.
 - Stone, John O. "Air pressure and cosmogenic isotope production." Journal of Geophysical Research: Solid Earth 105.B10 (2000): 23753-23759.
 - Titterington, D. Michael, and Alex N. Halliday. "On the fitting of parallel isochrons and the method of maximum likelihood." Chemical Geology 26.3-4 (1979): 183-195.
 - Vermeesch, Pieter. "IsoplotR: A free and open toolbox for geochronology." Geoscience Frontiers 9.5 (2018): 1479-1493.
 - William III, E. Odom. Dating the Cenozoic incision history of the Tennessee and Shenandoah Rivers with cosmogenic nuclides and 40Ar/39Ar in manganese oxides. Diss. Purdue University Graduate School, 2020.
 - York, Derek. "Least squares fitting of a straight line with correlated errors." Earth and planetary science letters 5 (1968): 320-324.
 - York, Derek, et al. "Unified equations for the slope, intercept, and standard errors of the best straight line." American journal of physics 72.3 (2004): 367-375.
