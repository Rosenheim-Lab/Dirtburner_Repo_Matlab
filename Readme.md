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

All files have help/doc text appended on to the front end of each m-file. Please use this text in each file to sort out what the input functions must be (all are .txt).

Devon testing