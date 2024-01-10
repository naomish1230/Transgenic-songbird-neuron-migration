# GitHub Repository (Transgenic-songbird-neuron-migration) for the following publication:

> In Vivo Imaging in Transgenic Songbirds Reveals Superdiffusive Neuron Migration in the Adult Brain (2024)
Naomi R. Shvedov, Sina Analoui, Theresia Dafalias, Timothy J. Gardner, Benjamin B. Scott


1. `/Transgenic-songbird-neuron-migration/Fig4-TrackVisualization/`
    * Code to plot the spatial positions and relative displacement of cell tracks in 3D, color-coded by time. 

2. `/Transgenic-songbird-neuron-migration/Fig5-CollectiveMigration/`
    * Data and code to analyze and plot pairwise correlations in speed and heading direction of all cells from multiple birds (`CollectiveMigrationFigureCorrelations.m`). 
    * Also includes a script to analyze and plot the entropy in position and  headings across the population of one bird (`CollectiveMigrationFigureEntropy.m`).

3. `/Transgenic-songbird-neuron-migration/Fig6-SuperdiffusiveMigration/`
    * All necessary code to recreate Figure 6. `Compute_MSDfunction.m` is called in `PopulationSuperdiffusivityPlots.m` script, which fits the diffusion model.

4. `/Transgenic-songbird-neuron-migration/Fig7-Simulation/`
    * This generative model is meant to simulate the migration of newborn neurons in HVC of juvenile songbirds as a stepwise process.
    * No additional file downloads are required, simulation can be downloaded and run as is.
    * Running this file yields a 3D representation of HVC that is intersected by the VZ with the birth and end points of simulated cells plotted in green and red, respectively.
    * Users can define how long the simulation is run, how many cells are being simulated, and if they want to plot the traces of each cell's migration path (See: User-Defined Parameters for Simulation Run)

5. `/Transgenic-songbird-neuron-migration/Functions/`
    * All functions that are called in above code.

6. `/Transgenic-songbird-neuron-migration/ProcessedData/`
    * Contains all data from experiments, both compiled (`combined_ALL`) and separated for each bird (`/IndividualBirdData/`). 
    * Data was tracked manually from registered and denoised 3D+t in vivo time-lapses taken with a two-photon microscope for >3 hours.
    * Data was pre-processed to remove cells with tracks that are too short, or are mis-clicks by the manual tracker.
    * Cells that have a tortuosity above 7 or a displacement below 5 were filtered out, as those were most likely tracked stationary cells.
    * `dynamic_stats.mat` files contain information for every cell position at every timepoint, including individual step distances and the real times of observation.
    * `neuroblast_stats.mat` files have the statistics for every cell track individually (e.g. total speed, tortuosity, total displacement, etc.).


## The data columns are labeled with abbreviations. The key is as follows: 

**FOR DYNAMIC STATS**

* N --> neuron ID
* I --> Intensity (a.u.)
* X,Y,Z --> Location in X,Y,Z respectively, in microns
* Slice --> the specific depth in the imaging Z-stack
* F --> the "frame" or timepoint of observation
* T --> extracted timestamp from metadata
* Q --> data was split up into 4 quadrants for easier tracking
* R --> region. 1 = HVC. 0 = hyperpallium/nidopallium
* S --> session. number of sessions varies across birds.
* MaxX,Y,Z,F --> X,Y,Z,F dimensions of whole quadrant
* Bird_ID --> shorthand bird identifier
* Age --> age of bird
* Sex --> 1 = male, 0 = female
* fixedQ --> 1 = quadrants were restitched after separation, XYZ positions are accurate relative to each other for all cells
* H --> How many hours have elapsed since first timepoint of that cell
* pathdiffs, pathXdiffs,pathYdiffs,pathZdiffs --> step size (Euclidean 3D distance between current and last position and X,Y,Z components, respectively)
* XY/XZangles --> step-wise heading direction in radians
* delt_XY/XZangles --> change in angle between current and previous step
* thresh... --> same as above except after a step distance threshold of 5 is applied 
* instvel --> instantaneous velocity, path distance over time since last step
* currentdisp/X/Y/Z --> cumulative displacement from starting position in 3D and in X,Y,Z components


**FOR NEUROBLAST STATS**

* pathtot (and pathX/Y/Z) --> total distance traveled by cell, both in 3D and in X/Y/Z components
* meanZpos --> average depth of cell
* duration --> length of observation in hours
* ...
* (see above)
* ...
* disptot (and dispX/Y/Z) --> total displacement of cell relative to start, both in 3D and in X/Y/Z components
* torts --> cell tortuosity
* speed --> total distance/total time
* dispXY/XZangle --> overall heading of displacement vector
* std/avg_threshdeltXY/XZ* --> average or St. Dev of thresholded (see above) changes in angle of all steps of cell
    * If cell travels in a straight line and does not deviate, the average thresh. deltXY/XZ will be 0, and St. Dev would be small