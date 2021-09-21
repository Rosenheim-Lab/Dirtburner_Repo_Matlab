# Rosenheim Laboratory Matlab Files for Processing Ramped PyrOx (RPO) Data
These files comprise the suite of m-files that allow you to:
1. visualize raw RPO data, 
2. plot raw RPO data from different runs together,
3. correct RPO data for both modern and dead blank contamination, 
4. plot and visualize RPO isotope data in a standardized way.

The files and their dependencies can be summarized as follows: 
* Visualize raw RPO data
   - DecompDB, DecompDB2, and DecompDB3. These functions ask for a file path and a vector of temperatures at which the samples were defined (empty vectors are accepted). The 3 different functions treat data from different stages in the development of the data acquisition system of the Dirtburner, and the help text at the beginning of each file (accessed through `help(func)` or `doc(func)`) tells the user which run series can be used with each function. No dependencies.
* Plot raw RPO data from different runs together
   - the DecompDB* series of functions outputs a normalized dataset of temperature and pCO_{2}. These series can fill different output namespaces for different runs and can be manually plotted together.
   - the DecompDB* series of functions can also operate through a directory of runs in a for loop, but remember that each iteration of any of the functions will output 3 graphs. Consider surpressing graphing output when using in this way.
* Correct RPO data for blank contamination
   - LiveBlankDeterm and DeadBlankDeterm both take in the blank data from previous runs of the RPO system and calculates an average and a standard deviation in blank mass for live blank (Fm = 1), and dead blank (Fm = 0, respectively). These functions need inputs of a filepath as well as a confidence interval. There is a graphical output (3-4 figures) to visualize the time evolution and analyst dependency on blank contamination. 
* Correct RPO data for blank contamination
   - BlankCorrect14C and BlankCorrect14_comps both calculate new fractions modern and uncertainties for input data. The functions call for a matrix of data (manually input) as well as live blank mass and uncertainty and dead blank mass and uncertainty. The output is a matrix with new fractions modern and uncertainties. 
* PlotThermo, PlotThermo_D14C, and PlotThermo_Fm all plot the thermographs with uneven bar circles to depict the decomposition reaction in light of isotope values. A plot is generated for visualizaton. The function calles two files - data (use _alt suffixed data for a better plot) and results (containing fractions modern, uncertainty, and stable isotope values as well as the temperature intervals of sampling). These functions differ only in how they visualize with uneven bars (positive and negative). Dependencies are with fracmat.m and unevenbarcircles.m. 
* Lierhaugen0-30 is a script that performs Gaussian decomposition of the reaction for sample Lierhaugen0_30_AT.txt and isotope results Lierhaugen0_30_DeminResults. The dependencies of this script (dfdp.m and nlleasqu.m) are also in the repository so the script should operate. To apply this script to a different set of run data and isotope results, one would have to:
  1. change the file names for the run data and the isotope results,
  2. determine number of Gaussians to model,
  3. change the the initial guess (pin) Gaussian values, and 
     * Each Gaussian contains three variables - height, center, width
     * Thus you need n*3 variables 
     * Width is tricky because it represents one standard deviation, or 67% of the area under the Gaussian. This will plot wider than you intended, so pick about 50% of the width you    suspect.
     * If you enter a number of pin values that is not a multiple of 3, you will get an error.
  4. change the function form to match your number of Gaussians - IMPORTANT!
  5. test the results until the script produces a stable convergence
Generally, this routine is very forgiving of "bad" initial guesses. It will converge with a good model. If it takes more than about 10 iterations, you should evaluate whether you are using the correct number Gaussians - your eye is really good at picking them out, don't overthink it!
This script can be used on dummy isotope values if none are available (shape runs or fast ramps), or you can excerpt the part of the script that does the Gaussian decomposition without running the isotope quantifications of Gaussian curves. 

All files have help/doc text appended on to the front end of each m-file. Please use this text in each file to sort out what the input functions must be (all are .txt).
